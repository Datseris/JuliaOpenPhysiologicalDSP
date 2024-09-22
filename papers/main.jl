using DrWatson
@quickactivate "JuliaOpenPhysiologicalDSP"
using CairoMakie
using ComplexityMeasures
using CSV
using DSP
# using DataFrames
using MAT
# using TopoPlots

# themeing
ENV["BGCOLOR"] = "white"
using MakieThemeing

# %% Load
"""
References for creating this parser come from here:

https://sccn.ucsd.edu/wiki/Makoto%27s_useful_EEGLAB_code#How_to_load_EEGLAB_.set_and_.fdt_files_without_using_EEGLAB_.2805.2F09.2F2020_updated.29
"""
function fdt_parser(fdt_file, mat_file)
    open(fdt_file) do io
        read!(
            io,
            Array{Float32}(undef, (Int(mat_file["nbchan"]), Int(mat_file["pnts"]))),
        )
    end
end
sub_1_ses_1_path = datadir("exp_raw", "sub-01", "ses-1", "eeg")
sub_2_ses_1_path = datadir("exp_raw", "sub-02", "ses-1", "eeg")
sub_1_mat_file = matread(joinpath(sub_1_ses_1_path, "sub-01_ses-1_task-eyesclosed_eeg.set"))
sub_2_mat_file = matread(joinpath(sub_2_ses_1_path, "sub-02_ses-1_task-eyesclosed_eeg.set"))
sub_1_data = fdt_parser(joinpath(sub_1_ses_1_path, "sub-01_ses-1_task-eyesclosed_eeg.fdt"), sub_1_mat_file)
sub_2_data = fdt_parser(joinpath(sub_2_ses_1_path, "sub-02_ses-1_task-eyesclosed_eeg.fdt"), sub_2_mat_file)

# actal timeseries
x = sub_1_data[1, :]
y = sub_2_data[1, :]
fs = 500.0 # sampling frequency

# subsample them
x = x[1:10:end]
y = y[1:10:end]
fs = 5000.0

t = range(0.002, 300; step = 0.02) # can we obtain this from the data?

# %% example
fig, axs = axesgrid(2, 1;
    size = (1200, 800), xlabels = "Time (sec)", sharex = true,
    ylabels = "Voltage (μV)"
)

lines!(axs[1], t, x)
lines!(axs[2], t, y; color = Cycled(2))

axs[1].title = "Subject 1 Raw Signal (Freq: 500 Hz)"
axs[2].title = "Subject 2 Raw Signal (Freq: 500 Hz)"

fig

# %% power bands

function decompose_into_bands(x, fs = 500.0)
    delta_range = [1, 4]
    theta_range = [4, 8]
    alpha_range = [8, 12]
    beta_range = [12, 30]
    gamma_rage = [30, 200]
    ranges = [delta_range, theta_range, alpha_range, beta_range, gamma_rage]
    signals = map(ranges) do r
        Ws = (r[1] - 0.5, r[2] + .5)./(fs/2)
        Wp = (r[1], r[2])./(fs/2)
        Rp = 3  # Passband ripple
        Rs = 40 # Stopband attenuation
        n, Wn = buttord(Wp, Ws, Rp, Rs)
        btr_band_filt = digitalfilter(Bandpass(r[1], r[2]; fs = fs), Butterworth(n))
        signal = filt(btr_band_filt, x)
    end
    pushfirst!(signals, x)
    return signals
end

signalsx = decompose_into_bands(x)
signalsy = decompose_into_bands(y)

# %%

fig, axs = axesgrid(length(signalsx), 2; sharex = true, sharey = true,
size = (700, 800),
    titles = ["Subject $i" for i in 1:2],
    ylabels = ["Voltage, $(b)" for b in ["all", "δ", "θ", "α", "β", "γ"]]
)

for (j, ts) in enumerate((signalsx, signalsy))
    for (i, s) in enumerate(ts)
        lines!(axs[i, j], t, s)
        xlims!(axs[i, j], 10, 50)
        ylims!(axs[i, j], nothing, nothing)
    end
end

fig


# %%
using DrWatson
@quickactivate "JuliaOpenPhysiologicalDSP"

using ComplexityMeasures
using DelayEmbeddings
using NeuroAnalyzer

function decompose_eeg(set_file)

function decompose_eeg(set_file; channel = 2, time = 20)
    eeg = import_set(set_file);
    NeuroAnalyzer.filter!(eeg, fprototype=:fir, ftype=:hp, cutoff=0.5, ch = eeg.locs[!, 1])
    s, _, _ = bpsplit(eeg, ch = eeg.locs[!, 1])

    full = eeg.data[2, 1:time*sr(eeg), 1]
    delta = s[1, channel, 1:time*sr(eeg), 1]
    theta = s[2, channel, 1:time*sr(eeg), 1]
    alpha = s[3, channel, 1:time*sr(eeg), 1]
    beta = s[6, channel, 1:time*sr(eeg), 1]
    gamma = s[9, channel, 1:time*sr(eeg), 1]

end

sub_1_ses_1 = datadir("exp_raw", "sub-01", "ses-1", "eeg", "sub-01_ses-1_task-eyesopen_eeg.set");
sub_2_ses_1 = datadir("exp_raw", "sub-02", "ses-1", "eeg", "sub-02_ses-1_task-eyesopen_eeg.set");

signal_sub_1, signal_sub_2 = decompose_eeg.([sub_1_ses_1, sub_2_ses_1]);

measures = (
    (x, τ) -> entropy_normalized(OrdinalPatterns(; m = 3, τ), x),
    (x, τ) -> entropy_normalized(WaveletOverlap(), x),
    (x, τ) -> entropy_normalized(PowerSpectrum(), x),
    (x, τ) -> complexity_normalized(SampleEntropy(x; m = 2, τ), x),
)

mnames = ["Perm-3", "Wavelet", "Spectral", "Sample", "ReverseD"]

fig, axs = axesgrid(1, 2; sharex = true, sharey = true,
    size = (700, 200), ylabels = "complexity",
    titles = ["Subject $i" for i in 1:2],
)

bands_idxs = [1, 2, 3, 4, 5]
band_names = ["δ", "θ", "α", "β", "γ"]

for (j, ts) in enumerate((signal_sub_1, signal_sub_2))
    ax = axs[j]
    for (k, m) in enumerate(measures)
        for (n, i) in enumerate(bands_idxs)
            s = ts[i]
            τ = estimate_delay(s, "mi_min", 1:100)
            spacing = 1/(length(bands_idxs)+1)
            q = m(s, τ)
            loc = k*1 + (n - 1)*spacing
            barplot!(ax, loc, q; width = spacing, color = Cycled(n))
            text!(ax, loc, q; text = band_names[n], align = (:center, :bottom))
        end
    end
    ax.xticks = (1:length(mnames), mnames)
    ylims!(ax, 0, 1.2)
    hideydecorations!(ax; label = false)
    ax.xgridvisible = false
end

display(fig)

# wsave(papersdir("main.png"), fig)
