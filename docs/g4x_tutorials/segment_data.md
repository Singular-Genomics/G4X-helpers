<br>

# Segmenting G4X (Cellpose)

Here, we lay out a basic approach to resegment your G4X data. We strongly encourage you to start this process by looking at the existing segmentation using our [G4X Viewer](../g4x_viewer/index.md). This will allow you to identify where it fails to capture your desired cell morphology. 

These instructions are designed to help you take the .jp2 images described in the [G4X data output](../g4x_data/g4x_output.md) documentation and run your own segmentation using Cellpose 2. This is typically an iterative process.

<br>
---

### Basic segmentation approach

For this tutorial, we use the most basic Cellpose arguments. It is strongly encouraged that you look into the [Cellpose documentation](https://cellpose.readthedocs.io/en/latest/) for more advanced options that may suit your needs.

#### 1. Import libraries
```python
import numpy as np
import glymur
from cellpose import models
from cellpose.models import CellposeModel
from skimage.exposure import rescale_intensity
```
<br>

#### 2. Choose which type of segmentation you want to run

=== "Nuclear Only"

    #### 3. Load the .jp2 image using glymur
    ```python
    # --- Load JP2 Images ---
    nuclear = glymur.Jp2k("path/to/nuclear.jp2")[:]
    ```

    #### 4. Normalize your nuclear image to the range 0-1
    ```python
    # --- Normalize nuclear image ---
    nuclear_norm = rescale_intensity(nuclear, out_range=(0, 1))
    ```

    #### 5. Run nuclear segmentation
    ```python
    # --- Run nuclear model ---
    model = models.CellposeModel(gpu=True, model_type="nuclear")
    masks, flows, styles = model.eval(nuclear_norm_subset, channels=[0],
                                        normalize=False)
    ```

=== "Nuclear +1 marker"

    #### 3. Load the .jp2 images using glymur
    ```python
    # --- Load JP2 Images ---
    nuclear = glymur.Jp2k("path/to/nuclear.jp2")[:]
    marker1 = glymur.Jp2k("path/to/marker1.jp2")[:]
    ```

    #### 4. Normalize each image channel to the range 0-1
    ```python
    # --- Normalize each channel separately ---
    nuclear_norm = rescale_intensity(nuclear, out_range=(0, 1))
    marker1_norm = rescale_intensity(marker1, out_range=(0, 1))
    ```
    
    #### 5. Stack your channels for Cellpose input
    ```python
    # --- Stack channels:  marker1_norm = channel 0, nuclear = channel 1 ---
    input_img = np.stack([marker1_norm, nuclear_norm], axis=2)
    ```

    #### 6. Run Cyto3 segmentation
    ```python
    # --- Run Cyto3 ---
    model = models.Cellpose(gpu=True, model_type="cyto3")
    masks, flows, styles, diams = model.eval(input_img, channels=[1, 0], 
                                                normalize=False)
    ```

=== "Nuclear +2 markers"
    
    #### 3. Load the .jp2 images using glymur
    ```python
    # --- Load JP2 Images ---
    nuclear = glymur.Jp2k("path/to/nuclear.jp2")[:]
    marker1 = glymur.Jp2k("path/to/marker1.jp2")[:]
    marker2 = glymur.Jp2k("path/to/marker2.jp2")[:]
    ```

    #### 4. Normalize each image channel to the range 0-1
    ```python
    # --- Normalize each channel separately ---
    nuclear_norm = rescale_intensity(nuclear, out_range=(0, 1))
    marker1_norm = rescale_intensity(marker1, out_range=(0, 1))
    marker2_norm = rescale_intensity(marker2, out_range=(0, 1))
    ```
    
    #### 5. Perform averaging of markers to create a combined membrane channel
    === "Simple Average"
        ```python
        # combine membrane markers
        membrane_combo = (marker1_norm + marker2_norm) / 2  # simple average
        ```
    === "Weighted Average"
        ```python
        # weighted average (sum of weights should be 1)
        marker1_weight = 0.35
        marker2_weight = 0.65
        membrane_combo = ((marker1_weight * marker1_norm) + 
                            (marker2_weight * marker2_norm))
        ```
    
    #### 6. Stack your channels for Cellpose input
    ```python
    # --- Stack channels:  membrane = channel 0, nuclear = channel 1 ---
    input_img = np.stack([membrane_combo, nuclear_norm], axis=2)
    ```
    #### 7. Run Cyto3 segmentation
    ```python
    # --- Run Cyto3 ---
    model = models.Cellpose(gpu=True, model_type="cyto3")
    masks, flows, styles, diams = model.eval(input_img, channels=[1, 0], 
                                                normalize=False)
    ```
---

<br>

### Iterating on your segmentation

After you've generated your segmentation masks, you sould take time to visualize them and determine if they are meeting your needs. This is an iterative process that may take several rounds of tuning. 

To save time, you may want to consider using some basic python plotting as well as subsetting your data to speed up iteration. Below are examples of each.

---
<br>

#### Visualize your segmentation

Here is a snippet of code to visualize your segmentation masks overlaid on a single fluorescent image. Notably, the variable this uses are:

* `marker_image`: a single fluorescent image (e.g. nuclear or membrane) to show the cell outlines against (same dimensions as `masks`)
* `masks`: the generated segmentation masks from Cellpose

```python
import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(8,8))

# Show the fluorescent image in blue with a bit of transparency
vmax = np.percentile(marker_image, 99) # remove hot pixels
plt.imshow(
    marker_image, 
    cmap="Blues", vmin=0, vmax=vmax, alpha=0.8, origin="lower"
)

# Draw one contour around each labeled cell
# using half-integer levels to outline integer labels cleanly
levels = np.arange(0.5, masks.max() + 0.5, 1.0)
plt.contour(masks, levels=levels, colors="gray", 
                linewidths=0.7, origin="lower") 

plt.axis("off")
plt.tight_layout()
plt.show()

```

<br>

#### Subset your images

Since this process can take quite a while on our stitched images, you may want to try to hone your approach first on a smaller, cropped image. For example, after loading your .jp2 image(s), you can subset them with something like this, running your segmentation on that block only:

```python

sz = 1000 # max size of your subset block

# get the image size
x = marker_norm.shape[0]
y = marker_norm.shape[1]

# take the middle 1000 x 1000 pixel block
marker_subset = marker_norm[int((y/2)-int(sz/2)):int((y/2)+int(sz/2)),
                            int((x/2)-int(sz/2)):int((x/2)+int(sz/2))]
```

### Applying your segmentation

Once you are happy with the outputs, you are ready to save them and run our G4X-helpers tool to regenerate your single cell data using your new segmentation masks.

---

#### Save your segmentation masks

```python
# ---- Save generated masks ----
np.save("/path/to/output/folder/resegmented_masks.npy", masks)
np.savez_compressed("/path/to/output/folder/resegmented_masks.npz", 
                    nuclei=masks, nuclei_exp=masks)
```

#### Regenerate your G4X outputs

After you have a segmentation that you are happy with, use out resegmentation tool to regenerate it for each tissue block. For installation tips, see our guide on [installing G4X-helpers](../g4x_helpers/installation/index.md). For advanced usage tips, see our [resegment documentation](../g4x_helpers/features/resegment.md).

```bash
resegment 
    --run_base /path/to/data/dir/<block_id>
    --segmentation_mask /path/to/new/masks.npy 
    --out_dir /path/to/output/dir/
```

<br>

--8<-- "_partials/end_cap.md"
