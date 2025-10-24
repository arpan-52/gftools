# Complete Guide: Processing uGMRT Full-Stokes Data

## Setup & Installation

### Step 1: Install Conda
Follow: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

### Step 2: Create Isolated Environment
```bash
conda create -n time_series python=3.12
conda activate time_series
```
This keeps your system clean and prevents conflicts with other projects.

### Step 3: Install Required Libraries
```bash
pip install your numpy matplotlib scipy
```

**What you're installing:**
- `your` - Read/write filterbank files (astronomy-specific I/O library)
- `numpy` - Numerical computing
- `matplotlib` - Plotting and visualization  
- `scipy` - Scientific functions

---

## Converting Raw GMRT Data to Filterbank Format

### Input Files
Download these files from GMRT (all needed):
```
J1932+1059_bm1_pa_550_200_16_22may2025.raw          (main data file)
J1932+1059_bm1_pa_550_200_16_22may2025.raw.hdr      (GMRT internal, not bothered)
J1932+1059_bm1_pa_550_200_16_22may2025.raw.bhdr     (main header)
J1932+1059_bm1_pa_550_200_16_22may2025.raw.ahdr     (start time and date in IST)
J1932+1059_bm1_pa_550_200_16_22may2025.raw.h4k      (GMRT internal, not bothered)
J1932+1059_bm1_pa_550_200_16_22may25_043407.log     (observation log)
```

### Convert to Stokes Parameters
The raw data is in a special format. Convert to readable filterbank files using the custom code, this package will be published soon:
```bash
python3 gmrt_convert.py J1932+1059_bm1_pa_550_200_16_22may2025.raw -t stokes
```

The `-t stokes` flag outputs 4 separate files for each polarization component.

### Output Files Created
```
J1932p1059_I.fil     (Total Intensity)
J1932p1059_Q.fil     (Linear component 1)
J1932p1059_U.fil     (Linear component 2)
J1932p1059_V.fil     (Circular component)
```

Each file is a 2D time-frequency array with the same dimensions but different polarization data.

---

## Understanding Your Data

### Check File Properties
Open Jupyter notebook or Python:

```python
from your import Your
import numpy as np
import matplotlib.pyplot as plt

# Load one file to inspect
your_obj = Your("J1932p1059_I.fil")

# Display header info
print(your_obj.your_header)
```

### Key Header Parameters Explained

```python
print(f"fch1: {your_obj.your_header.fch1}")           # Output: 550.0 (MHz)
print(f"foff: {your_obj.your_header.foff}")           # Output: 0.09765625 (MHz)
print(f"nchans: {your_obj.your_header.nchans}")       # Output: 2048 (channels)
print(f"tsamp: {your_obj.your_header.tsamp}")         # Output: 0.00016384 (seconds)
print(f"nspectra: {your_obj.your_header.nspectra}")   # Output: 1835008 (samples)
```

| Parameter | Value | Meaning |
|-----------|-------|---------|
| **fch1** | 550.0 MHz | Highest frequency channel |
| **foff** | 0.0977 MHz | Frequency spacing between channels (goes 550, 550.098, 550.196...) |
| **nchans** | 2048 | Total frequency channels |
| **tsamp** | 0.000164 s | Time between samples |
| **nspectra** | 1,835,008 | Total time samples (~300 seconds) |

**Data shape:** (1,835,008 time samples, 2048 frequency channels)

---

## Loading & Inspecting Data

### Load All Stokes Components
```python
I_obj = Your("J1932p1059_I.fil")
Q_obj = Your("J1932p1059_Q.fil")
U_obj = Your("J1932p1059_U.fil")
V_obj = Your("J1932p1059_V.fil")

# Get data arrays
I = I_obj.get_data(nstart=0, nsamp=I_obj.your_header.nspectra)
Q = Q_obj.get_data(nstart=0, nsamp=Q_obj.your_header.nspectra)
U = U_obj.get_data(nstart=0, nsamp=U_obj.your_header.nspectra)
V = V_obj.get_data(nstart=0, nsamp=V_obj.your_header.nspectra)

print(f"Data shape: {I.shape}")  # (1835008, 2048)
```

### If Running Out of Memory
Load smaller chunks:
```python
nsamp = 100000  # Load first 100k samples instead of all
I = I_obj.get_data(nstart=0, nsamp=nsamp)
Q = Q_obj.get_data(nstart=0, nsamp=nsamp)
U = U_obj.get_data(nstart=0, nsamp=nsamp)
V = V_obj.get_data(nstart=0, nsamp=nsamp)
```

---

## Dedispersion: Why and How

### The Problem
Radio signals travel through space containing free electrons. These electrons delay lower frequencies more than higher frequencies. When you observe across a wide frequency range, the signal gets "smeared out" in time - it looks dispersed.

### The Solution
Dedispersion corrects this by shifting each frequency channel back by the right amount so all frequencies arrive aligned.

### The Formula
```
Time Delay = 4.15 ms × DM × [(ν₁/GHz)⁻² - (ν₂/GHz)⁻²]
```
Where:
- **DM** = Dispersion Measure (pc/cm³) - tells us how many electrons between us and source
- **ν** = Frequency in GHz

### Implement Dedispersion Function
```python
def dedisperse(data, your_obj, dm):
    """
    Dedisperse data at a given dispersion measure.
    
    Args:
        data: time-frequency array (time_samples, frequency_channels)
        your_obj: Your object containing header info
        dm: Dispersion measure in pc/cm³
    
    Returns:
        Dedispersed array (same shape as input)
    """
    
    # Step 1: Create frequency array in GHz
    freq_ghz = (your_obj.your_header.fch1 + 
                np.arange(your_obj.your_header.nchans) * 
                your_obj.your_header.foff) / 1000
    
    # Step 2: Calculate time delay for each frequency (in milliseconds)
    delay_ms = 4.15 * dm * (1/freq_ghz**2 - 1/freq_ghz[0]**2)
    
    # Step 3: Convert delay from ms to sample indices
    delays = (delay_ms / (your_obj.your_header.tsamp * 1000)).astype(int)
    
    # Step 4: Apply shifts to each frequency channel
    dedispersed = np.zeros_like(data)
    for i in range(data.shape[1]):
        dedispersed[:, i] = np.roll(data[:, i], -delays[i])
    
    return dedispersed
```

### Apply Dedispersion
```python
dm = 560.0  # Dispersion measure (from literature or search)

I_d = dedisperse(I, I_obj, dm)
Q_d = dedisperse(Q, Q_obj, dm)
U_d = dedisperse(U, U_obj, dm)
V_d = dedisperse(V, V_obj, dm)

print(f"Dedispersed data shape: {I_d.shape}")
```

---

## Computing Polarization Parameters

From the 4 Stokes parameters (I, Q, U, V), we derive:

```python
# Linear polarization (Q and U components combined)
L = np.sqrt(Q_d**2 + U_d**2)

# Circular polarization (V component)
C = np.abs(V_d)

# Polarization fraction as percentage
# (add 1e-10 to avoid dividing by zero)
L_fraction = (L / (I_d + 1e-10)) * 100
C_fraction = (C / (I_d + 1e-10)) * 100
```

**Interpretation:**
- **L_fraction = 0%** → Signal is unpolarized
- **L_fraction = 100%** → Signal is fully linearly polarized
- **C_fraction = 0%** → No circular polarization
- **C_fraction = 100%** → Fully circularly polarized

---

## Visualization

### Plot 1: All Stokes Parameters
```python
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# I: Total intensity
axes[0, 0].imshow(I_d.T, aspect='auto', cmap='viridis', origin='lower')
axes[0, 0].set_title('Stokes I (Total Intensity)', fontweight='bold')
axes[0, 0].set_ylabel('Frequency Channel')

# Q: Linear polarization component 1
axes[0, 1].imshow(Q_d.T, aspect='auto', cmap='RdBu', origin='lower')
axes[0, 1].set_title('Stokes Q', fontweight='bold')

# U: Linear polarization component 2
axes[1, 0].imshow(U_d.T, aspect='auto', cmap='RdBu', origin='lower')
axes[1, 0].set_title('Stokes U', fontweight='bold')
axes[1, 0].set_xlabel('Time Sample')
axes[1, 0].set_ylabel('Frequency Channel')

# V: Circular polarization
axes[1, 1].imshow(V_d.T, aspect='auto', cmap='RdBu', origin='lower')
axes[1, 1].set_title('Stokes V (Circular)', fontweight='bold')
axes[1, 1].set_xlabel('Time Sample')

fig.suptitle(f'Stokes Parameters (Dedispersed at DM={dm})', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('01_stokes_parameters.png', dpi=300)
plt.show()
```

**What to look for:**
- **I**: Should show a clear pulse (bright spot in time-frequency plane)
- **Q, U**: May show polarization structure; can be positive or negative
- **V**: Shows circular polarization; usually weaker than Q, U

### Plot 2: Polarization Fractions
```python
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Linear polarization percentage
im1 = axes[0].imshow(L_fraction.T, aspect='auto', cmap='hot', 
                      origin='lower', vmin=0, vmax=100)
axes[0].set_title('Linear Polarization Fraction (%)', fontweight='bold')
axes[0].set_xlabel('Time Sample')
axes[0].set_ylabel('Frequency Channel')
cbar1 = plt.colorbar(im1, ax=axes[0], label='Percentage')

# Circular polarization percentage
im2 = axes[1].imshow(C_fraction.T, aspect='auto', cmap='plasma', 
                      origin='lower', vmin=0, vmax=100)
axes[1].set_title('Circular Polarization Fraction (%)', fontweight='bold')
axes[1].set_xlabel('Time Sample')
cbar2 = plt.colorbar(im2, ax=axes[1], label='Percentage')

fig.suptitle(f'Polarization Analysis (DM={dm})', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('02_polarization_fractions.png', dpi=300)
plt.show()
```

**Interpretation:**
- **Hot colors** (red) = High polarization
- **Cool colors** (dark) = Low/unpolarized emission

---

## Complete Analysis Script

Save and run this for full processing:

```python
from your import Your
import numpy as np
import matplotlib.pyplot as plt

print("=" * 50)
print("uGMRT Full-Stokes Analysis")
print("=" * 50)

# LOAD DATA
print("\n[1/5] Loading data...")
I_obj = Your("J1932p1059_I.fil")
Q_obj = Your("J1932p1059_Q.fil")
U_obj = Your("J1932p1059_U.fil")
V_obj = Your("J1932p1059_V.fil")

I = I_obj.get_data(nstart=0, nsamp=I_obj.your_header.nspectra)
Q = Q_obj.get_data(nstart=0, nsamp=Q_obj.your_header.nspectra)
U = U_obj.get_data(nstart=0, nsamp=U_obj.your_header.nspectra)
V = V_obj.get_data(nstart=0, nsamp=V_obj.your_header.nspectra)
print(f"✓ Loaded data shape: {I.shape}")

# DEDISPERSION FUNCTION
print("\n[2/5] Setting up dedispersion...")
def dedisperse(data, your_obj, dm):
    freq_ghz = (your_obj.your_header.fch1 + 
                np.arange(your_obj.your_header.nchans) * 
                your_obj.your_header.foff) / 1000
    delay_ms = 4.15 * dm * (1/freq_ghz**2 - 1/freq_ghz[0]**2)
    delays = (delay_ms / (your_obj.your_header.tsamp * 1000)).astype(int)
    dedispersed = np.zeros_like(data)
    for i in range(data.shape[1]):
        dedispersed[:, i] = np.roll(data[:, i], -delays[i])
    return dedispersed

# DEDISPERSE
print("\n[3/5] Dedispersing...")
dm = 560.0
I_d = dedisperse(I, I_obj, dm)
Q_d = dedisperse(Q, Q_obj, dm)
U_d = dedisperse(U, U_obj, dm)
V_d = dedisperse(V, V_obj, dm)
print(f"✓ Dedispersion complete at DM={dm} pc/cm³")

# CALCULATE POLARIZATION
print("\n[4/5] Computing polarization...")
L = np.sqrt(Q_d**2 + U_d**2)
C = np.abs(V_d)
L_frac = (L / (I_d + 1e-10)) * 100
C_frac = (C / (I_d + 1e-10)) * 100
print("✓ Polarization computed")

# PLOTTING
print("\n[5/5] Generating plots...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes[0, 0].imshow(I_d.T, aspect='auto', cmap='viridis', origin='lower')
axes[0, 0].set_title('Stokes I')
axes[0, 1].imshow(Q_d.T, aspect='auto', cmap='RdBu', origin='lower')
axes[0, 1].set_title('Stokes Q')
axes[1, 0].imshow(U_d.T, aspect='auto', cmap='RdBu', origin='lower')
axes[1, 0].set_title('Stokes U')
axes[1, 1].imshow(V_d.T, aspect='auto', cmap='RdBu', origin='lower')
axes[1, 1].set_title('Stokes V')
fig.suptitle(f'Stokes Parameters (DM={dm})')
plt.tight_layout()
plt.savefig('stokes.png', dpi=300)

fig2, axes2 = plt.subplots(1, 2, figsize=(14, 5))
im1 = axes2[0].imshow(L_frac.T, aspect='auto', cmap='hot', origin='lower', vmin=0, vmax=100)
axes2[0].set_title('Linear Polarization (%)')
plt.colorbar(im1, ax=axes2[0])
im2 = axes2[1].imshow(C_frac.T, aspect='auto', cmap='plasma', origin='lower', vmin=0, vmax=100)
axes2[1].set_title('Circular Polarization (%)')
plt.colorbar(im2, ax=axes2[1])
plt.tight_layout()
plt.savefig('polarization.png', dpi=300)

plt.show()
print("✓ Plots saved: stokes.png, polarization.png")
print("\n" + "=" * 50)
print("Analysis complete!")
print("=" * 50)
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Memory error | Reduce `nsamp` parameter when loading |
| Slow processing | Use smaller data chunk |
| NaN values in output | Check input data for zeros/NaNs |
| Black plots | Adjust colormap or use `vmin`/`vmax` |
