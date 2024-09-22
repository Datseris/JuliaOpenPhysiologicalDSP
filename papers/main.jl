using DrWatson
@quickactivate "JuliaOpenPhysiologicalDSP"

using CairoMakie
using ComplexityMeasures
using DelayEmbeddings
using NeuroAnalyzer

ENV["BGCOLOR"] = "white"
using MakieThemeing

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

    return [full, delta, theta, alpha, beta, gamma]
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

bands_idxs = [1, 2, 3, 4, 5, 6]
band_names = ["F", "δ", "θ", "α", "β", "γ"]

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

wsave(papersdir("main.png"), fig)
