#!/usr/bin/env python3
"""Exploratory CK18-associated fluorescence image maps: Eu versus Cy3.

NOT a calibrated expression map. Does not assert matching fields or comparable
photon intensities. Uses fixed per-image percentiles only for visualization.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, label
import tifffile

CHANNEL_INDEX = {'red': 0, 'green': 1, 'blue': 2}
EXTS = {'.tif', '.tiff', '.png', '.jpg', '.jpeg'}


def read_image(path: str | Path):
    p = Path(path)
    if p.suffix.lower() in {'.tif', '.tiff'}:
        data = tifffile.imread(p)
    else:
        data = np.asarray(Image.open(p))
    if data.ndim == 3 and data.shape[0] in (3, 4) and data.shape[-1] not in (3, 4):
        data = np.moveaxis(data, 0, -1)
    if data.ndim not in (2, 3) or (data.ndim == 3 and data.shape[-1] not in (3, 4)):
        raise ValueError(f'Unsupported image shape: {data.shape}')
    return data


def extract_channel(data, channel):
    if data.ndim == 2:
        return data.astype(np.float64)
    if channel == 'gray':
        raise ValueError('gray is only valid for single-channel images; choose red/green/blue')
    return data[..., CHANNEL_INDEX[channel]].astype(np.float64)


def display_scale(arr, low=1, high=99.5):
    a = np.asarray(arr, dtype=float)
    finite = a[np.isfinite(a)]
    if finite.size == 0:
        return np.zeros(a.shape, dtype=float)
    vmin, vmax = np.percentile(finite, [low, high])
    return np.clip((a-vmin) / max(float(vmax-vmin), 1e-12), 0, 1)


def original_preview(data):
    if data.ndim == 2:
        return display_scale(data)
    return np.stack([display_scale(data[..., n]) for n in range(3)], axis=-1)


def classify(path):
    parts = [p.lower() for p in path.parts]
    has_ck18 = any('ck18' in p for p in parts)
    if not has_ck18:
        return None
    modality = None
    for p in parts:
        if 'cy3' in p:
            modality = 'Cy3'
        if 'eu' in p and any(t in p for t in ('-eu','_eu','eu-','eu_')):
            modality = 'Eu'
    if modality is None:
        return None
    return modality


def scan_images(root, out, max_preview=36):
    root, out = Path(root), Path(out)
    if not root.is_dir():
        raise FileNotFoundError(root)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for p in sorted(root.rglob('*')):
        if not p.is_file() or p.suffix.lower() not in EXTS:
            continue
        modality = classify(p.relative_to(root))
        if not modality:
            continue
        try:
            img = read_image(p)
            desc = f'{img.shape} {img.dtype}'
            err = ''
        except Exception as e:
            img = None
            desc = ''
            err = str(e)[:160]
        rows.append({'modality': modality, 'relative_path': str(p.relative_to(root)),
                     'shape_dtype': desc, 'read_error': err,
                     'suggested_channel': 'red' if modality == 'Eu' else 'green',
                     'reviewed_same_marker': 'no', 'reviewed_signal_channel': 'no'})
    with (out / 'ck18_image_inventory.csv').open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]) if rows else ['relative_path'])
        w.writeheader(); w.writerows(rows)
    good = [r for r in rows if not r['read_error']]
    selected = []
    if good:
        # Interleave source directories to increase preview diversity.
        buckets = {}
        for item in good:
            key = (item['modality'], str(Path(item['relative_path']).parent))
            buckets.setdefault(key, []).append(item)
        while buckets and len(selected) < max_preview:
            for key in list(buckets):
                if buckets[key] and len(selected) < max_preview:
                    selected.append(buckets[key].pop(0))
                if not buckets[key]:
                    del buckets[key]
        ncol = 4
        nrow = int(np.ceil(len(selected)/ncol))
        fig, axarr = plt.subplots(nrow, ncol, figsize=(ncol*3.5, nrow*3.2), squeeze=False)
        for ax in axarr.flat:
            ax.axis('off')
        for ax, entry in zip(axarr.flat, selected):
            img = read_image(root / entry['relative_path'])
            ax.imshow(original_preview(img), cmap='gray' if img.ndim == 2 else None)
            ax.set_title(entry['modality'] + '  ' + str(Path(entry['relative_path']).parent), fontsize=8)
            ax.set_xlabel(Path(entry['relative_path']).name, fontsize=7)
            ax.axis('on'); ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle('CK18 candidate image inventory — CHECK modality, staining and channels', fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.98])
        fig.savefig(out / 'ck18_contact_sheet.png', dpi=130, bbox_inches='tight')
        plt.close(fig)
    print(f'Found {len(rows)} CK18 candidate images ({len(good)} readable).')
    print('Inventory:', out / 'ck18_image_inventory.csv')
    if good:
        print('Preview:', out / 'ck18_contact_sheet.png')
    return rows


def spatial_map(channel_data, sigma=12.0, mad_multiplier=3.0, border=0.05, min_area=6):
    """Local intensity contrast for display images, not a CK18 expression estimate."""
    arr = np.asarray(channel_data, dtype=np.float64)
    if arr.ndim != 2 or min(arr.shape) < 16:
        raise ValueError('Expected 2D image with both dimensions >=16')
    if not np.all(np.isfinite(arr)):
        raise ValueError('Nonfinite pixel values')
    if not 0 <= border < 0.40:
        raise ValueError('border must be >=0 and <0.4')
    h, w = arr.shape
    dy, dx = int(round(h*border)), int(round(w*border))
    valid = np.zeros((h,w), bool)
    valid[dy:h-dy if dy else h, dx:w-dx if dx else w] = True
    if valid.sum() < 20:
        raise ValueError('Insufficient valid pixels')
    reference = gaussian_filter(arr, sigma=sigma, mode='reflect')
    residual = arr-reference
    values = residual[valid]
    center = float(np.median(values))
    mad = float(np.median(np.abs(values-center)))
    robust_sigma = max(1.4826*mad, 1e-8)
    threshold = center + mad_multiplier * robust_sigma
    initial = (residual > threshold) & valid
    # Remove isolated tiny components; no segmentation semantics.
    comp, count = label(initial)
    if count:
        sizes = np.bincount(comp.ravel())
        mask = initial & (sizes[comp] >= max(1, min_area))
    else:
        mask = initial
    positive = np.maximum(residual, 0.)
    contrast_disp = positive / max(float(np.percentile(positive[valid], 99.5)), 1e-8)
    contrast_disp = np.clip(contrast_disp, 0, 1)
    return {'raw': arr, 'reference': reference, 'residual': residual,
            'relative_map': contrast_disp, 'valid': valid, 'mask': mask,
            'threshold': float(threshold), 'robust_sigma': float(robust_sigma),
            'fraction': float(mask[valid].mean()),
            'residual_p95': float(np.percentile(residual[valid], 95)),
            'residual_median': float(np.median(residual[valid])),
            'selected_pixel_count': int(mask.sum()),
            'valid_pixel_count': int(valid.sum())}


def save_single(img, computed, image_name, modality, channel, outdir):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 4, figsize=(16, 4.5))
    ax[0].imshow(original_preview(img), cmap='gray' if img.ndim == 2 else None)
    ax[0].set_title('Original exported image')
    ax[1].imshow(display_scale(computed['raw']), cmap='gray')
    ax[1].set_title(f'Extracted {channel} channel')
    m = ax[2].imshow(computed['relative_map'], cmap='magma', vmin=0, vmax=1)
    ax[2].set_title('Within-image local contrast')
    fig.colorbar(m, ax=ax[2], fraction=0.043, label='Display-normalized')
    ax[3].imshow(display_scale(computed['raw']), cmap='gray')
    overlay = np.zeros((*computed['mask'].shape, 4), dtype=float)
    overlay[computed['mask']] = (1.0, 0.15, 0.12, 0.80)
    ax[3].imshow(overlay)
    ax[3].set_title(f'High-contrast candidates: {computed["fraction"]*100:.2f}%')
    for a in ax: a.axis('off')
    fig.suptitle(f'{modality} — {image_name} | CK18-associated image intensity proxy ONLY', fontsize=12)
    fig.text(.5, .015, 'Not an absolute CK18 expression map, positive-cell mask, or calibrated photon image.',
             ha='center', color='dimgray', fontsize=9)
    fig.tight_layout(rect=[0, .06, 1, .93])
    output=outdir / f'{modality.lower()}_ck18_signal_map.png'
    fig.savefig(output, dpi=170, bbox_inches='tight')
    plt.close(fig)
    return output


def compare(cy3, eu, cy3_channel, eu_channel, out, sigma, border, mad_multiplier):
    out=Path(out); out.mkdir(parents=True, exist_ok=True)
    metrics=[]; maps={}
    for modality, filename, channel in [('Cy3', cy3, cy3_channel), ('Eu', eu, eu_channel)]:
        p=Path(filename)
        if not p.exists(): raise FileNotFoundError(p)
        data=read_image(p)
        a=extract_channel(data, channel)
        q=spatial_map(a, sigma=sigma, border=border, mad_multiplier=mad_multiplier)
        save_single(data, q, p.name, modality, channel, out)
        maps[modality]=(data,q,p)
        metrics.append({'modality':modality, 'file':str(p), 'channel':channel,
                        'shape':str(data.shape), 'dtype':str(data.dtype),
                        'high_local_contrast_percent':round(q['fraction']*100, 4),
                        'local_residual_p95_display_units':round(q['residual_p95'],4),
                        'threshold_display_units':round(q['threshold'],4),
                        'selected_pixels':q['selected_pixel_count'],
                        'valid_pixels':q['valid_pixel_count']})
    with (out/'ck18_image_metrics.csv').open('w', newline='', encoding='utf-8-sig') as f:
        writer=csv.DictWriter(f,fieldnames=list(metrics[0]))
        writer.writeheader();writer.writerows(metrics)
    fig, axarr=plt.subplots(2, 3, figsize=(12, 8))
    for row, modality in enumerate(['Cy3', 'Eu']):
        data, q, p = maps[modality]
        axarr[row,0].imshow(original_preview(data), cmap='gray' if data.ndim==2 else None)
        axarr[row,0].set_title(f'{modality}: original ({p.name})', fontsize=10)
        axarr[row,1].imshow(q['relative_map'], cmap='magma', vmin=0,vmax=1)
        axarr[row,1].set_title(f'{modality}: within-image local contrast', fontsize=10)
        axarr[row,2].imshow(display_scale(q['raw']),cmap='gray')
        overlay=np.zeros((*q['mask'].shape,4),dtype=float)
        overlay[q['mask']] = [1,.12,.12,.82]
        axarr[row,2].imshow(overlay)
        axarr[row,2].set_title(f'{modality}: candidate high contrast ({q["fraction"]*100:.2f}%)', fontsize=9)
        for a in axarr[row]: a.axis('off')
    fig.suptitle('CK18-associated signal maps | Eu and Cy3 are independently normalized', fontsize=13)
    fig.text(.5,.015,'Different fields/acquisition/color encodings: no absolute intensity, expression, or fold-change comparison.',
             ha='center', fontsize=9, color='dimgray')
    fig.tight_layout(rect=[0,.04,1,.95])
    fig.savefig(out/'ck18_eu_cy3_side_by_side.png', dpi=170, bbox_inches='tight')
    plt.close(fig)
    metadata={'interpretation':'exploratory image-intensity proxies only, not CK18-positive cell classification',
              'side_by_side_not_registered': True, 'physical_quantitative_eu_cy3_comparison': False,
              'cy3_channel':cy3_channel, 'eu_channel':eu_channel,
              'gaussian_sigma_px':sigma,'border_fraction_per_side':border,'mad_multiplier':mad_multiplier,
              'signal_background_note':'Background is a smoothed local image, not an experimentally verified negative control'}
    (out/'analysis_notes.json').write_text(json.dumps(metadata,indent=2,ensure_ascii=False),encoding='utf-8')
    print('Completed:',out)
    for r in metrics:
        print(f"{r['modality']} | {r['channel']} | exploratory high-local-contrast area: {r['high_local_contrast_percent']:.2f}%")
    print('Figure:',out/'ck18_eu_cy3_side_by_side.png')
    print('CAUTION: no calibrated Eu/Cy3 fold-change or biological CK18 expression claims.')
    return metrics


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    subs=parser.add_subparsers(dest='command',required=True)
    s=subs.add_parser('scan',help='Inventory and contact sheet for CK18-labelled folders')
    s.add_argument('--root',required=True)
    s.add_argument('--out',required=True)
    s.add_argument('--max-preview',type=int,default=36)
    c=subs.add_parser('compare',help='Create independent image maps for one Cy3 and one Eu CK18 file')
    c.add_argument('--cy3-image',required=True)
    c.add_argument('--eu-image',required=True)
    c.add_argument('--cy3-channel',choices=['red','green','blue','gray'],default='green')
    c.add_argument('--eu-channel',choices=['red','green','blue','gray'],default='red')
    c.add_argument('--sigma',type=float,default=12)
    c.add_argument('--border',type=float,default=0.05)
    c.add_argument('--mad-multiplier',type=float,default=3.0)
    c.add_argument('--out',required=True)
    a=parser.parse_args()
    if a.command=='scan': scan_images(a.root,a.out,a.max_preview)
    else: compare(a.cy3_image,a.eu_image,a.cy3_channel,a.eu_channel,a.out,a.sigma,a.border,a.mad_multiplier)


if __name__=='__main__':
    main()
