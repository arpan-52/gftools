# SOP for using full-Stokes uGMRT data

## First, we need to setup a nice working environment with all the softwares

We will be using a set of custom made tools and also public softwares to test things out.

### Installing Softwares

First install conda:
Follow the instructions on https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

When successfully installed, it can be called as conda from terminal.

Now creating an environment so that we secure us and also others.

In the terminal, use:
```bash
conda create -n time_series python=3.12
```

-n sets the name of the environment and python=3.12 decides the python version to use.

Upon completion we want to activate the environment.
```bash
conda activate time_series
```

Now we want to install some software, in the terminal use:
```bash
pip3 install your
```

This will install github.com/thepetabyteproject/your/

This is our main I/O code for filterbank files.

Now, the fun part, we will be using uGMRT raw files which is internal and not readable by any other software, we will use our own custom code, shared with the email itself. This is a custom code and will be published as a package till now, we will use the python code itself.

## Download the uGMRT data from the link

The data will come as this raw data and all the headers, we need all these files to run our code and it is mostly provided by the GMRT itself.

```
J1932+1059_bm1_pa_550_200_16_22may2025.raw.hdr
J1932+1059_bm1_pa_550_200_16_22may2025.raw.bhdr
J1932+1059_bm1_pa_550_200_16_22may2025.raw.ahdr
J1932+1059_bm1_pa_550_200_16_22may25_043407.log
J1932+1059_bm1_pa_550_200_16_22may2025.raw.h4k
J1932+1059_bm1_pa_550_200_16_22may2025.raw
```

## To convert the raw file to filterbank

We are going to use gftools package, the code is:
```bash
python3 gmrt_convert.py J1932+1059_bm1_pa_550_200_16_22may2025.raw -t stokes
```

-t stokes will basically output I Q U V as 4 separate Stokes filterbanks. The raw data contain a sandwiched format from all the correlations.

Now you can also do `python3 gmrt_convert.py -h` to show all the options.

This step will create:
```
J1932p1059_I.fil
J1932p1059_Q.fil
J1932p1059_U.fil
J1932p1059_V.fil
```

## Now, we will plot and dedisperse them to see how they look like

We either can use Jupyter notebook or Python console. Here are some example codes.

First we import the libraries in your Jupyter notebook or Python console:

```python
from your import Your
import numpy as np
import matplotlib.pyplot as plt
```

Load one filterbank file:

```python
your_obj = Your("J1932p1059_Q.fil")
data = your_obj.get_data(nstart=0, nsamp=your_obj.your_header.nspectra)
```

Check the header parameters:

```python
print(f"fch1: {your_obj.your_header.fch1}")
print(f"foff: {your_obj.your_header.foff}")
print(f"nchans: {your_obj.your_header.nchans}")
print(f"tsamp: {your_obj.your_header.tsamp}")
```

Output:
```
fch1: 550.0
foff: 0.09765625
nchans: 2048
tsamp: 0.00016384
```

**What these parameters mean:**
- `fch1`: First frequency channel is 550 MHz
- `foff`: Frequency offset between channels is 0.098 MHz (spacing)
- `nchans`: Total number of frequency channels (2048)
- `tsamp`: Time between samples in seconds (sampling interval)

Now we have all the necessary params.

## So we will now define a dedispersion routine

The dispersion measure (DM) corrects for the delay caused by free electrons in space. Different frequencies arrive at different times. We need to shift each frequency channel back by the correct amount so they align.

The formula is:
```
time_delay = 4.15 ms × DM × [(ν₁/GHz)⁻² - (ν₂/GHz)⁻²]
```

Define the dedispersion function in your notebook:

```python
def dedisperse(data, your_obj, dm):
    """
    Dedisperse filterbank data at a given DM
    
    Parameters:
    data - time-frequency array (nsamples, nchans)
    your_obj - Your object with header info
    dm - dispersion measure (pc/cm³)
    """
    # Convert frequencies to GHz
    freq_ghz = (your_obj.your_header.fch1 + 
                np.arange(your_obj.your_header.nchans) * 
                your_obj.your_header.foff) / 1000
    
    # Calculate delay in ms for each frequency
    delay_ms = 4.15 * dm * (1/freq_ghz**2 - 1/freq_ghz[0]**2)
    
    # Convert to sample indices
    delays = (delay_ms / (your_obj.your_header.tsamp * 1000)).astype(int)
    
    # Shift each frequency channel by its delay
    dedispersed = np.zeros_like(data)
    for i in range(data.shape[1]):
        dedispersed[:, i] = np.roll(data[:, i], -delays[i])
    
    return dedispersed
```

## Load all Stokes parameters

```python
I_obj = Your("J1932p1059_I.fil")
Q_obj = Your("J1932p1059_Q.fil")
U_obj = Your("J1932p1059_U.fil")
V_obj = Your("J1932p1059_V.fil")

I = I_obj.get_data(nstart=0, nsamp=I_obj.your_header.nspectra)
Q = Q_obj.get_data(nstart=0, nsamp=Q_obj.your_header.nspectra)
U = U_obj.get_data(nstart=0, nsamp=U_obj.your_header.nspectra)
V = V_obj.get_data(nstart=0, nsamp=V_obj.your_header.nspectra)
```

## Dedisperse all Stokes parameters

```python
dm = 560.0  # Dispersion measure

I_d = dedisperse(I, I_obj, dm)
Q_d = dedisperse(Q, Q_obj, dm)
U_d = dedisperse(U, U_obj, dm)
V_d = dedisperse(V, V_obj, dm)
```

## Calculate polarization parameters

From Stokes parameters we can derive:
- Linear polarization: `L = √(Q² + U²)`
- Circular polarization: `C = |V|`
- Polarization fraction: percentage of polarized signal

```python
# Linear polarization
L = np.sqrt(Q_d**2 + U_d**2)

# Circular polarization
C = np.abs(V_d)

# As percentages of total intensity
L_frac = (L / (I_d + 1e-10)) * 100
C_frac = (C / (I_d + 1e-10)) * 100
```

## Now plot I Q U V filterbanks

```python
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

axes[0, 0].imshow(I_d.T, aspect='auto', cmap='viridis', origin='lower')
axes[0, 0].set_title('Stokes I (Intensity)')
axes[0, 0].set_ylabel('Frequency Channels')

axes[0, 1].imshow(Q_d.T, aspect='auto', cmap='RdBu', origin='lower')
axes[0, 1].set_title('Stokes Q')

axes[1, 0].imshow(U_d.T, aspect='auto', cmap='RdBu', origin='lower')
axes[1, 0].set_title('Stokes U')
axes[1, 0].set_xlabel('Time Samples')
axes[1, 0].set_ylabel('Frequency Channels')

axes[1, 1].imshow(V_d.T, aspect='auto', cmap='RdBu', origin='lower')
axes[1, 1].set_title('Stokes V (Circular)')
axes[1, 1].set_xlabel('Time Samples')

plt.suptitle(f'Stokes Parameters - Dedispersed at DM={dm}', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('stokes_parameters.png', dpi=300)
plt.show()
```

## Then plot linear and circular polarization as percentage

```python
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

im1 = axes[0].imshow(L_frac.T, aspect='auto', cmap='hot', origin='lower', vmin=0, vmax=100)
axes[0].set_title('Linear Polarization Fraction (%)')
axes[0].set_xlabel('Time Samples')
axes[0].set_ylabel('Frequency Channels')
plt.colorbar(im1, ax=axes[0], label='Percentage (%)')

im2 = axes[1].imshow(C_frac.T, aspect='auto', cmap='plasma', origin='lower', vmin=0, vmax=100)
axes[1].set_title('Circular Polarization Fraction (%)')
axes[1].set_xlabel('Time Samples')
axes[1].set_ylabel('Frequency Channels')
plt.colorbar(im2, ax=axes[1], label='Percentage (%)')

plt.suptitle(f'Polarization Analysis - DM={dm}', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('polarization_fractions.png', dpi=300)
plt.show()
```

## Mention them if the code does not work, read only some samples

If you run into memory issues or the code is too slow, read only a subset of the data:

```python
nsamp = 100000  # Read only 100k samples instead of all
I = I_obj.get_data(nstart=0, nsamp=nsamp)
Q = Q_obj.get_data(nstart=0, nsamp=nsamp)
U = U_obj.get_data(nstart=0, nsamp=nsamp)
V = V_obj.get_data(nstart=0, nsamp=nsamp)
```

Then continue with dedispersion and plotting on this smaller dataset.

## Complete script

Save this as a Python file or run in Jupyter notebook:

```python
from your import Your
import numpy as np
import matplotlib.pyplot as plt

# Load data
I_obj = Your("J1932p1059_I.fil")
Q_obj = Your("J1932p1059_Q.fil")
U_obj = Your("J1932p1059_U.fil")
V_obj = Your("J1932p1059_V.fil")

I = I_obj.get_data(nstart=0, nsamp=I_obj.your_header.nspectra)
Q = Q_obj.get_data(nstart=0, nsamp=Q_obj.your_header.nspectra)
U = U_obj.get_data(nstart=0, nsamp=U_obj.your_header.nspectra)
V = V_obj.get_data(nstart=0, nsamp=V_obj.your_header.nspectra)

# Define dedispersion
def dedisperse(data, your_obj, dm):
    freq_ghz = (your_obj.your_header.fch1 + np.arange(your_obj.your_header.nchans) * your_obj.your_header.foff) / 1000
    delay_ms = 4.15 * dm * (1/freq_ghz**2 - 1/freq_ghz[0]**2)
    delays = (delay_ms / (your_obj.your_header.tsamp * 1000)).astype(int)
    dedispersed = np.zeros_like(data)
    for i in range(data.shape[1]):
        dedispersed[:, i] = np.roll(data[:, i], -delays[i])
    return dedispersed

# Dedisperse
dm = 560.0
I_d = dedisperse(I, I_obj, dm)
Q_d = dedisperse(Q, Q_obj, dm)
U_d = dedisperse(U, U_obj, dm)
V_d = dedisperse(V, V_obj, dm)

# Calculate polarization
L = np.sqrt(Q_d**2 + U_d**2)
C = np.abs(V_d)
L_frac = (L / (I_d + 1e-10)) * 100
C_frac = (C / (I_d + 1e-10)) * 100

# Plot Stokes
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes[0, 0].imshow(I_d.T, aspect='auto', cmap='viridis', origin='lower')
axes[0, 0].set_title('Stokes I')
axes[0, 1].imshow(Q_d.T, aspect='auto', cmap='RdBu', origin='lower')
axes[0, 1].set_title('Stokes Q')
axes[1, 0].imshow(U_d.T, aspect='auto', cmap='RdBu', origin='lower')
axes[1, 0].set_title('Stokes U')
axes[1, 1].imshow(V_d.T, aspect='auto', cmap='RdBu', origin='lower')
axes[1, 1].set_title('Stokes V')
plt.suptitle(f'Stokes Parameters - DM={dm}')
plt.tight_layout()
plt.savefig('stokes.png', dpi=300)
plt.show()

# Plot polarization fractions
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
im1 = axes[0].imshow(L_frac.T, aspect='auto', cmap='hot', origin='lower', vmin=0, vmax=100)
axes[0].set_title('Linear Polarization (%)')
plt.colorbar(im1, ax=axes[0])
im2 = axes[1].imshow(C_frac.T, aspect='auto', cmap='plasma', origin='lower', vmin=0, vmax=100)
axes[1].set_title('Circular Polarization (%)')
plt.colorbar(im2, ax=axes[1])
plt.tight_layout()
plt.savefig('polarization.png', dpi=300)
plt.show()
```
