#!/usr/bin/env python3
"""Chapter 2 Eu/Cy3 IHC contrast audit and explicit-ROI descriptive measurements.

This script NEVER assumes TIFF RGB pixels are calibrated emission intensities.
All ROI positions use normalized x,y image coordinates [0,1].
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import tifffile
from PIL import Image, ImageDraw, ImageOps

IMAGE_EXTS = {'.tif', '.tiff'}
XLS_EXTS = {'.xlsx', '.xls', '.xlsm'}
ROI_KEYS = ('signal_x0', 'signal_y0', 'signal_x1', 'signal_y1', 'bg_x0', 'bg_y0', 'bg_x1', 'bg_y1')
MANIFEST_FIELDS = ['image_id','cell_line','modality','replicate_id','image_path','channel',*ROI_KEYS,'roi_verified','intensity_provenance']


def read_image(path: Path) -> np.ndarray:
    a = tifffile.imread(path)
    if a.ndim == 2 or (a.ndim == 3 and a.shape[-1] in (3,4)):
        return a
    raise ValueError(f'Unsupported image shape {a.shape}; requires 2D, RGB, or RGBA')


def channel_plane(a: np.ndarray, channel: str) -> np.ndarray:
    if a.ndim == 2:
        return a.astype(np.float64)
    if channel not in ('red','green','blue'):
        raise ValueError('For RGB/RGBA TIFF select --channel red/green/blue')
    return a[..., ('red','green','blue').index(channel)].astype(np.float64)


def make_preview(a: np.ndarray, max_size=(390, 260)) -> Image.Image:
    # Preview-only individual min/max stretch: NEVER use preview values as quantification input.
    x = a[..., :3] if a.ndim == 3 else a
    x = np.asarray(x, dtype=np.float32)
    small = x[::max(1,x.shape[0]//300),::max(1,x.shape[1]//300)]
    if x.ndim == 3:
        lo,hi = np.percentile(small, [1,99.5],axis=(0,1))
        rgb = np.clip((x-lo)/(hi-lo+1e-5),0,1)
        im = Image.fromarray(np.uint8(rgb*255),'RGB')
    else:
        lo,hi = np.percentile(small,[1,99.5])
        gray = np.uint8(np.clip((x-lo)/(hi-lo+1e-5),0,1)*255)
        im = Image.fromarray(gray,'L').convert('RGB')
    im.thumbnail(max_size, Image.Resampling.LANCZOS)
    return im


def modality_and_group(rel: Path) -> tuple[str,str,str]:
    folder = str(rel.parent)
    low = folder.lower()
    mode = 'Cy3' if 'cy3' in low else ('Eu' if 'eu' in low else 'unknown')
    line = 'BT474' if re.search(r'(?:^|[/\\])474(?:[/\\]|$)',low) or 'bt474' in low else ('MCF7' if 'mcf7' in low else 'unverified')
    return mode,line,folder


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open('w',newline='',encoding='utf-8-sig') as f:
        w = csv.DictWriter(f,fieldnames=fields,extrasaction='ignore')
        w.writeheader();w.writerows(rows)


def scan(args) -> None:
    root = Path(args.root).expanduser().resolve()
    if not root.is_dir(): raise SystemExit(f'Data directory not found: {root}')
    out = Path(args.out).expanduser().resolve();out.mkdir(parents=True,exist_ok=True)
    files = sorted(p for p in root.rglob('*') if p.is_file() and p.suffix.lower() in (IMAGE_EXTS|XLS_EXTS) and out not in p.parents)
    rows=[];xls=[];representatives={};counts=Counter();manifest=[]
    for p in files:
        rel=p.relative_to(root)
        modality,line,group=modality_and_group(rel)
        if p.suffix.lower() in XLS_EXTS:
            xls.append(str(rel));continue
        row={'image_path':str(rel),'modality_inferred':modality,'cell_line_inferred':line,'group':group,'file_mb':round(p.stat().st_size/1048576,2)}
        try:
            with tifffile.TiffFile(p) as t:
                series=t.series[0]
                shp=series.shape
                row.update(shape=str(shp),dtype=str(series.dtype),height=shp[0],width=shp[1],channels=(shp[-1] if len(shp)==3 else 1))
            if len(representatives.get(group,[]))<args.previews_per_group:
                representatives.setdefault(group,[]).append(p)
        except Exception as e:
            row['read_error']=str(e)
        rows.append(row);counts[(modality,line)]+=1
    write_csv(out/'image_inventory.csv',rows,['image_path','modality_inferred','cell_line_inferred','group','file_mb','shape','dtype','height','width','channels','read_error'])
    (out/'spreadsheets_found.txt').write_text(('\n'.join(xls) if xls else 'No Excel files found under root')+'\n',encoding='utf-8')
    tiles=[]
    for group,paths in sorted(representatives.items()):
        for p in paths:
            if len(tiles)>=args.max_previews:break
            tile=Image.new('RGB',(430,330),'white')
            try:
                im=make_preview(read_image(p),(410,255))
                tile.paste(im,((430-im.width)//2,8))
                d=ImageDraw.Draw(tile)
                d.text((9,270),f'{group[:52]}',fill='black')
                d.text((9,290),f'{p.name[:55]}',fill='black')
                tiles.append(tile)
            except Exception as e:
                print(f'WARN: preview {p.name}: {e}')
        if len(tiles)>=args.max_previews:break
    if tiles:
        cols=3; sheet=Image.new('RGB',(cols*430,((len(tiles)+cols-1)//cols)*330),'#e8e8e8')
        for i,tile in enumerate(tiles): sheet.paste(tile,((i%cols)*430,(i//cols)*330))
        sheet.save(out/'contact_sheet.png')
    # Blank ROI manifest: avoid pretending these are matched views or verified biological replicates.
    for r in rows:
        if r['modality_inferred'] not in ('Eu','Cy3') or r['cell_line_inferred'] not in ('BT474','MCF7'):continue
        if r.get('read_error'):continue
        manifest.append({'image_id':Path(r['image_path']).stem+'_'+str(len(manifest)),
                         'cell_line':r['cell_line_inferred'],'modality':r['modality_inferred'],
                         'replicate_id':'UNKNOWN','image_path':r['image_path'],'channel':'REVIEW',
                         'roi_verified':'no','intensity_provenance':'unknown'})
    write_csv(out/'roi_manifest_template.csv',manifest,MANIFEST_FIELDS)
    lines=[f'Image files: {len(rows)}', f'Excel files: {len(xls)}',f'Preview images: {len(tiles)}','', 'Images by inferred group:']
    lines += [f'  {key[0]} / {key[1]}: {val}' for key,val in sorted(counts.items())]
    lines += ['', 'All cell-type and modality labels are inferred from filenames; verify experimentally.',
              'TIFF stored RGB values may be color-rendered/quantized; not calibrated photons.',
              'Spreadsheet filenames are listed, but workbook contents have NOT been interpreted.']
    (out/'summary.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('\n'.join(lines));print(f'\nSaved: {out}')


def inspect(args) -> None:
    p=Path(args.image).expanduser()
    a=read_image(p)
    c=channel_plane(a,args.channel)
    fig,ax=plt.subplots(figsize=(9,6))
    lo,hi=np.percentile(c,[1,99.5]); ax.imshow(c,vmin=lo,vmax=max(hi,lo+1e-4),cmap='gray',extent=(0,1,1,0))
    ax.set_xticks(np.arange(0,1.01,0.1));ax.set_yticks(np.arange(0,1.01,0.1))
    ax.grid(color='cyan',alpha=.4,linewidth=.6)
    ax.set_title(f'{p.name} | {args.channel} | {a.shape} — DISPLAY ONLY')
    ax.set_xlabel('normalized x: left 0 → right 1');ax.set_ylabel('normalized y: top 0 → bottom 1')
    fig.tight_layout()
    out=Path(args.out).expanduser();out.parent.mkdir(parents=True,exist_ok=True);fig.savefig(out,dpi=160);plt.close(fig)
    print(f'Saved: {out}')
    print('Use grid to select two separate, nonoverlapping ROIs; edit manifest CSV normalized coordinates.')


def roi_mask(shape: tuple[int,int], xyxy: tuple[float,float,float,float]) -> np.ndarray:
    h,w=shape; x0,y0,x1,y1=(float(v) for v in xyxy)
    if not (0<=x0<x1<=1 and 0<=y0<y1<=1):raise ValueError(f'Invalid ROI [0,1] xyxy: {xyxy}')
    xa=int(np.floor(x0*w));xb=int(np.ceil(x1*w))
    ya=int(np.floor(y0*h));yb=int(np.ceil(y1*h))
    mask=np.zeros((h,w),dtype=bool);mask[ya:yb,xa:xb]=True
    return mask


def quantify(plane: np.ndarray, srect, brect) -> dict:
    s=roi_mask(plane.shape,srect);b=roi_mask(plane.shape,brect)
    if np.any(s&b):raise ValueError('Signal and background ROIs overlap')
    a=plane[s];d=plane[b]
    if len(a)<9 or len(d)<9:raise ValueError('ROIs require at least 9 pixels each')
    sm=float(a.mean());bm=float(d.mean());bs=float(d.std(ddof=1))
    return {'signal_mean':sm,'background_mean':bm,'background_sd':bs,
            'delta_mean':sm-bm,'SBR':sm/bm if bm>0 else float('nan'),
            'CNR':(sm-bm)/bs if bs>0 else float('nan'),
            'n_signal_pixels':len(a),'n_background_pixels':len(d)}


def measure(args) -> None:
    root=Path(args.root).expanduser().resolve();out=Path(args.out).expanduser().resolve()
    out.mkdir(parents=True,exist_ok=True);(out/'qc').mkdir(exist_ok=True)
    with Path(args.manifest).expanduser().open(encoding='utf-8-sig',newline='') as f: records=list(csv.DictReader(f))
    results=[];skips=[]
    for i,r in enumerate(records):
        name=r.get('image_id',str(i))
        if str(r.get('roi_verified','')).strip().lower() not in ('yes','true','1'):
            continue
        try:
            if r.get('channel') not in ('red','green','blue'):
                raise ValueError('Select channel red/green/blue after channel verification')
            if r.get('modality') not in ('Eu','Cy3'):
                raise ValueError('modality must be Eu or Cy3')
            vals=[r[k] for k in ROI_KEYS]
            if any(x is None or str(x).strip()=='' for x in vals):raise ValueError('Missing ROI coordinates')
            rect_s=tuple(map(float,vals[:4]));rect_b=tuple(map(float,vals[4:]))
            path=(root/r['image_path']).resolve()
            if not path.is_relative_to(root):raise ValueError('Image outside root')
            a=read_image(path);c=channel_plane(a,r['channel'])
            metrics=quantify(c,rect_s,rect_b)
            res={k:r.get(k,'') for k in ('image_id','cell_line','modality','replicate_id','image_path','channel','intensity_provenance')}
            res.update(metrics);results.append(res)
            # Display-only histogram and boxes (measurement uses ORIGINAL TIFF numeric array)
            fig,ax=plt.subplots(figsize=(8,5))
            lo,hi=np.percentile(c,[1,99.5]);ax.imshow(c,vmin=lo,vmax=max(hi,lo+1e-4),cmap='gray',extent=(0,1,1,0))
            import matplotlib.patches as patches
            for rect,color,lab in ((rect_s,'lime','signal'),(rect_b,'yellow','background')):
                x0,y0,x1,y1=rect
                ax.add_patch(patches.Rectangle((x0,y0),x1-x0,y1-y0,fill=False,edgecolor=color,linewidth=2,label=lab))
            ax.legend(loc='lower right');ax.set_title(f'{name}: CNR={metrics["CNR"]:.3f}, SBR={metrics["SBR"]:.3f}\nExploratory; RGB values may be display-encoded')
            fig.tight_layout();fname=re.sub(r'[^A-Za-z0-9_.-]+','_',name)[:80]
            fig.savefig(out/'qc'/f'{fname}_roi_qc.png',dpi=160);plt.close(fig)
            print(f'OK  {name}: SBR={metrics["SBR"]:.3f} CNR={metrics["CNR"]:.3f}')
        except Exception as e:
            skips.append({'image_id':name,'error':str(e)});print(f'SKIP {name}: {e}')
    fields=['image_id','cell_line','modality','replicate_id','image_path','channel','intensity_provenance',
            'signal_mean','background_mean','background_sd','delta_mean','SBR','CNR','n_signal_pixels','n_background_pixels']
    write_csv(out/'roi_measurements.csv',results,fields)
    write_csv(out/'skipped.csv',skips,['image_id','error'])
    # Descriptive plot only, never unqualified Eu/Cy3 advantage claim.
    if results:
        labels=sorted(set((r['cell_line'],r['modality']) for r in results))
        fig,ax=plt.subplots(figsize=(max(7,1.45*len(labels)),4.8))
        rng=np.random.default_rng(0)
        for idx,(line,mode) in enumerate(labels):
            values=[r['CNR'] for r in results if r['cell_line']==line and r['modality']==mode and np.isfinite(r['CNR'])]
            if values:ax.scatter(idx+rng.uniform(-.07,.07,len(values)),values,label='per image' if idx==0 else None)
        ax.set_xticks(range(len(labels)),[f'{line}\n{mode}' for line,mode in labels])
        ax.set_ylabel('CNR (ROI mean difference / background SD)')
        ax.set_title('Exploratory within-image ROI contrast\nNOT evidence of Eu superiority without validated acquisition and independent samples')
        ax.axhline(0,color='gray',lw=.7);fig.tight_layout()
        fig.savefig(out/'contrast_exploratory.png',dpi=180);plt.close(fig)
    (out/'interpretation.txt').write_text(
        'Within-image ROI SBR and CNR. NOT published SNR unless definitions are confirmed.\n'
        'No pixel is a biological replicate. Compare modalities only with matched biological replicates,\n'
        'validated target/background ROIs, identical processing definitions, and documented image linearity/normalization.\n'
        'RGB/RGBA TIFF exports may be display-encoded and cannot establish calibrated Eu vs Cy3 advantage.\n',encoding='utf-8')
    print(f'Saved: {out}. Measured {len(results)} verified-ROI images; skipped {len(skips)} errors. Unverified rows ignored.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('scan',help='Audit TIFF and Excel filenames; preview representative TIFFs')
    s.add_argument('--root',required=True);s.add_argument('--out',required=True)
    s.add_argument('--previews-per-group',type=int,default=1);s.add_argument('--max-previews',type=int,default=16)
    s.set_defaults(func=scan)
    s=sub.add_parser('inspect',help='Show normalized coordinate grid for choosing ROI boxes')
    s.add_argument('--image',required=True);s.add_argument('--channel',choices=['red','green','blue'],required=True)
    s.add_argument('--out',required=True);s.set_defaults(func=inspect)
    s=sub.add_parser('measure',help='Quantify verified, nonoverlapping ROI boxes from manifest')
    s.add_argument('--root',required=True);s.add_argument('--manifest',required=True);s.add_argument('--out',required=True)
    s.set_defaults(func=measure)
    args=p.parse_args();args.func(args)

if __name__=='__main__':main()
