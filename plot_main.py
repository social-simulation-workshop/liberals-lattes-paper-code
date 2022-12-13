from plot import PlotLinesHandler
import numpy as np

if __name__ == "__main__":
    NETSIZE_INDEX = 5
    ABSCORR_INDEX = 10
    NETSIZE_INTERVAL = 500

    res = [[] for _ in range(10)]
    with open("results/outcome_figure3_on_sta.txt", "r") as f:
        for l in f.read().split("\n"):
            if l:
                l = list(map(float, l.split(" ")))
                res[int(l[NETSIZE_INDEX]/NETSIZE_INTERVAL)-1].append(l[ABSCORR_INDEX])
    res = np.array(res)
    res = res[:, ~np.isnan(res).any(axis=0)]
    res_mean = np.array([np.mean(r) for r in res])
    res_std = np.array([np.std(r) for r in res])
                
    fn_suffix = "N_500_5000_std"
    plot_handler = PlotLinesHandler(xlabel="Network Size",
                                    ylabel="Mean Pairwise Correlation Magnitude",
                                    title=None,
                                    fn="fig3",
                                    x_lim=[500, 5000], x_tick=[500, 5000, 500], 
                                    y_lim=[0.0, 1.0], y_tick=[0.0, 1.0, 0.1],
                                    use_ylim=True,
                                    figure_ratio=853/1090,
                                    figure_size=7.5)
    plot_handler.plot_line_errorbar(res_mean, res_std, data_log_v=NETSIZE_INTERVAL, linewidth=2, color="black")
    plot_handler.save_fig(fn_suffix=fn_suffix)