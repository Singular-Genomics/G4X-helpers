<br>

# Plotting Spatial Data on the G4X

One of the first things you'll want to do with your spatial data is to visualize your data on a spatial plot. [G4X outputs](../g4x_data/g4x_output.md) are designed to make this information easy to access and plot in a variety of ways. With the tables and data structures provided, you can create a variety of spatial visualizations using tools like [Seurat](https://satijalab.org/seurat/articles/install_v5.html) in R or [scanpy](https://scanpy.readthedocs.io/en/stable/index.html) in Python.

<br>

### Python Visualization
---

This brief tutorial uses the scanpy library and assume that you have a working installation of it on your local system. If you do not, please see the [scanpy installation guide](https://scanpy.readthedocs.io/en/stable/getting_started/installation.html) for instructions on how to install it.


#### 1. Import libraries

```python
# Import necessary libraries
import scanpy as sc
import matplotlib.pyplot as plt
```

#### 2. Load your data

```python
# Load your data
adata = sc.read_h5ad('path/to/your/data.h5ad')
```

#### 3. Extract spatial coordinates into a scanpy `obsm` layer for scanpy compatibility

```python
# Extract spatial coordinates into obsm layer
adata.obsm['X_spatial_cell'] = adata.obs[['cell_x', 'cell_y']].to_numpy()
```
!!!note
    In order for scanpy to extract this information for default plotting, you have to put `(x, y)` coordinate pairs for each cell into the `obsm` layer and name it something starting with `X_`. You can do this with any pair of coordinates.

#### 4. Plot the spatial coordinates

```python
# Plot the spatial coordinates
plt.figure(figsize=(10,8))
ax = sc.pl.scatter(adata, basis='spatial_cell', color='log1p_total_counts', show=False)

# Flip the y-axis to match image orientation
ax.invert_yaxis()
plt.show()
```

!!!note
    Here, we're plotting with the `color='log1p_total_counts'` to visualize the log total transcript counts per cell, but you can replace this with any other variable in your `adata.obs` or `adata.var` that is the correct dimensions (1xncells) to visualize different aspects of your data.

#### 5. Save the plot as a PNG file

```python
# Save the plot as PNG
plt.savefig("spatial_plot.png", dpi=300, bbox_inches='tight')
```

<br>

### R Visualization
---

This brief tutorial uses the ggplot2 library and assumes that you have a working installation of it on your local system. If you do not, please see the [ggplot2 installation guide](https://ggplot2.tidyverse.org) for instructions on how to install it.

#### 1. Load necessary libraries

```R
# Load necessary libraries
library(ggplot2)
``` 

#### 2. Load your data and create a Seurat object

```R
# Load your data and create a Seurat object
run_base = c('/path/to/g4x_output/<sample_id>')
metadata_path = file.path(run_base, "single_cell_data/cell_metadata.csv.gz")
meta <- read.csv(gzfile(metadata_path), row.names = 1)
```

#### 3. Plot spatial coordinates with ggplot2

```R
# Ensure coordinates are numeric
meta$cell_x <- as.numeric(meta$cell_x)
meta$cell_y <- as.numeric(meta$cell_y)

p <- ggplot(meta, aes(x = cell_x, y = cell_y, color=log1p_total_counts)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_y_reverse() +           # <-- invert Y
  coord_equal() +
  labs(x = "cell_x", y = "cell_y", color = "log1p(total_counts)") +
  theme_minimal()

print(p)
```
!!!note
    Here, we're plotting with the `color='log1p_total_counts'` to visualize the log total transcript counts per cell, but you can replace this with any other variable in your `meta` dataframe that is the correct dimensions (1xncells) to visualize different aspects of your data.


#### 4. Save the plot as a PNG file

```R
# Save as Plot as PNG
ggplot2::ggsave(
  filename = "cell_xy_log1p_total_counts.png",
  plot     = p,
  width    = 8, height = 8, dpi = 300, units = "in",
  bg       = "white"   
)
```
<br>

 --8<-- "_partials/end_cap.md"
