import os
import datetime
import pandas as pd
import numpy as np
import sys
import re
import random
from idmtools.entities import IAnalyzer
from idmtools.entities.simulation import Simulation
import manifest

## For plotting
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates


class InsetChartAnalyzer(IAnalyzer):

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, expt_name, sweep_variables=None, channels=None, working_dir=".", start_year=0):
        super(InsetChartAnalyzer, self).__init__(working_dir=working_dir, filenames=["output/InsetChart.json"])
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.inset_channels = channels or ['Statistical Population', 'New Clinical Cases', 'Blood Smear Parasite Prevalence',
                                           'Infectious Vectors']
        self.expt_name = expt_name
        self.start_year = start_year

    def map(self, data, simulation: Simulation):
        simdata = pd.DataFrame({x: data[self.filenames[0]]['Channels'][x]['Data'] for x in self.inset_channels})
        simdata['Time'] = simdata.index
        simdata['Day'] = simdata['Time'] % 365
        simdata['Year'] = simdata['Time'].apply(lambda x: int(x / 365) + self.start_year)
        simdata['date'] = simdata.apply(
            lambda x: datetime.date(int(x['Year']), 1, 1) + datetime.timedelta(int(x['Day']) - 1), axis=1)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
            elif sweep_var == 'Run_Number':
                simdata[sweep_var] = 0
        return simdata

    def reduce(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        if not os.path.exists(os.path.join(self.working_dir, self.expt_name)):
            os.mkdir(os.path.join(self.working_dir, self.expt_name))

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv(os.path.join(self.working_dir, self.expt_name, 'All_Age_InsetChart.csv'), index=False)


class MonthlyPfPRAnalyzer(IAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir='./', start_year=0,
                 burnin=None, filter_exists=False):

        super(MonthlyPfPRAnalyzer, self).__init__(working_dir=working_dir,
                                                   filenames=["output/MalariaSummaryReport_monthly.json"]
                                                   )

        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name
        self.start_year = start_year
        self.burnin = burnin
        self.filter_exists = filter_exists

    def filter(self, simulation: Simulation):
        if self.filter_exists:
            file = os.path.join(simulation.get_path(), self.filenames[0])
            return os.path.exists(file)
        else:
            return True

    def map(self, data, simulation: Simulation):
        adf = pd.DataFrame()
        fname = self.filenames[0]
        age_bins = data[self.filenames[0]]['Metadata']['Age Bins']

        for age in range(len(age_bins)):
            d = data[fname]['DataByTimeAndAgeBins']['PfPR by Age Bin'][:-1]
            pfpr = [x[age] for x in d]

            d = data[fname]['DataByTimeAndAgeBins']['Annual Clinical Incidence by Age Bin'][:-1]
            clinical_cases = [x[age] for x in d]

            d = data[fname]['DataByTimeAndAgeBins']['Annual Severe Incidence by Age Bin'][:-1]
            severe_cases = [x[age] for x in d]

            d = data[fname]['DataByTimeAndAgeBins']['New Infections by Age Bin'][:-1]
            new_infections = [x[age] for x in d]

            d = data[fname]['DataByTimeAndAgeBins']['Annual Mild Anemia by Age Bin'][:-1]
            mild_anemia = [x[age] for x in d]

            d = data[fname]['DataByTimeAndAgeBins']['Pf Gametocyte Prevalence by Age Bin'][:-1]
            gamete_prev = [x[age] for x in d]

            d = data[fname]['DataByTimeAndAgeBins']['Mean Log Parasite Density by Age Bin'][:-1]
            dens = [x[age] for x in d]

            d = data[fname]['DataByTimeAndAgeBins']['Average Population by Age Bin'][:-1]
            pop = [x[age] for x in d]

            simdata = pd.DataFrame({'month': range(1, len(pfpr)+1),
                                    'PfPR': pfpr,
                                    'Cases': clinical_cases,
                                    'Severe_cases': severe_cases,
                                    'New_infections': new_infections,
                                    'Anemia': mild_anemia,
                                    'Mean_density': dens,
                                    'Gametocyte_prevalence': gamete_prev,
                                    'Pop': pop})

            simdata['agebin'] = age_bins[age]

            adf = pd.concat([adf, simdata])

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                try:
                    adf[sweep_var] = simulation.tags[sweep_var]
                except:
                    adf[sweep_var] = '-'.join([str(x) for x in simulation.tags[sweep_var]])

        adf['Year'] = None
        adf['Year'] = np.floor((adf['month']-1)/12)+self.start_year

        return adf

    def reduce(self, all_data):

        selected = [data for sim, data in all_data.items()]
        print(len(selected))
        if len(selected) == 0:
            print("\nWarning: No data have been returned... Exiting...")
            return

        if not os.path.exists(os.path.join(self.working_dir, self.expt_name)):
            os.mkdir(os.path.join(self.working_dir, self.expt_name))

        print(f'\nSaving outputs to: {os.path.join(self.working_dir, self.expt_name)}')

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv((os.path.join(self.working_dir, self.expt_name, 'PfPR_ClinicalIncidence_monthly.csv')),
                   index=False)


if __name__ == "__main__":
    # Lake Victoria seed sweep. Run:  python analyzer_W2.py
    # (analyzer_spatial.py must be in this folder)
    from idmtools.analysis.analyze_manager import AnalyzeManager
    from idmtools.core import ItemType
    from idmtools.core.platform_factory import Platform
    from analyzer_spatial import SpatialAnalyzer
    import glob

    EXP_ID = "63808b9e-03a6-4158-a1d2-f40d8cc47f67"     # your Lake Victoria sweep
    NAME = "lake_sweep_test"
    START_YEAR = 2010

    wdir = os.path.join(manifest.job_directory, "my_outputs", NAME)
    os.makedirs(wdir, exist_ok=True)
    sweep = ["Run_Number"]

    # ---- look at what the simulation actually wrote -----------------------------
    spatial_channels = ["Population", "PCR_Parasite_Prevalence", "New_Clinical_Cases"]
    have_spatial = False
    exp_dirs = [d for d in glob.glob(os.path.join(manifest.job_directory, "**", EXP_ID), recursive=True) if os.path.isdir(d)]
    print("Experiment folder(s):", exp_dirs)
    for ed in exp_dirs[:1]:
        outs = glob.glob(os.path.join(ed, "*", "output"))
        if outs:
            files = sorted(os.listdir(outs[0]))
            print("Output files of first simulation:\n ", "\n  ".join(files))
            have_spatial = all(os.path.exists(os.path.join(outs[0], f"SpatialReportMalariaFiltered_all_ages_{c}.bin")) for c in spatial_channels)
    print("Spatial files present:", have_spatial)
    channels = ["Statistical Population", "True Prevalence", "New Clinical Cases",
                "Infectious Vectors", "Rainfall", "Air Temperature"]

    with Platform("SLURM_LOCAL", job_directory=manifest.job_directory) as platform:
        analyzers = [
            InsetChartAnalyzer(expt_name=NAME, channels=channels, sweep_variables=sweep,
                               start_year=START_YEAR, working_dir=wdir),
            MonthlyPfPRAnalyzer(expt_name=NAME, sweep_variables=sweep,
                                start_year=START_YEAR, working_dir=wdir),
        ]
        if have_spatial:
            analyzers.append(SpatialAnalyzer(dir_name=NAME, f_base="SpatialReportMalariaFiltered", f_suffix="all_ages",
                            exp_id=EXP_ID, spatial_channels=spatial_channels,
                            sweep_variables=sweep, working_dir=wdir))
        else:
            print("\nNo spatial files in this experiment: skipping SpatialAnalyzer.")
        AnalyzeManager(configuration={}, ids=[(EXP_ID, ItemType.EXPERIMENT)],
                       analyzers=analyzers, partial_analyze_ok=True).analyze()

    if have_spatial:
        sp = pd.read_csv(os.path.join(wdir, "SpatialReportMalariaFiltered_all_ages.csv"))
        print("\nSpatial rows:", len(sp), " nodes:", sorted(sp["node"].unique()))
        print(sp.groupby("node")[["Population", "PCR_Parasite_Prevalence", "New_Clinical_Cases"]].mean().round(3))

        # ---- how much do the seeds differ? ---------------------------------------
        last = sp[sp["time"] >= sp["time"].max() - 180]            # last 6 months of the year
        seed_tab = last.groupby(["node", "Run_Number"])[["PCR_Parasite_Prevalence", "New_Clinical_Cases"]].mean().round(4)
        piv = seed_tab["PCR_Parasite_Prevalence"].unstack("Run_Number")
        piv["mean"] = piv.mean(axis=1)
        piv["sd"] = piv.drop(columns="mean").std(axis=1)
        print("\nPCR prevalence by county and seed (mean of last 6 months):")
        print(piv.round(3))
        seed_tab.to_csv(os.path.join(wdir, "seed_noise_by_county.csv"))
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, axs = plt.subplots(2, 4, figsize=(14, 6), sharex=True, sharey=True)
            for ax, (node, g) in zip(axs.ravel(), sp.groupby("node")):
                for run, gg in g.groupby("Run_Number"):
                    ax.plot(gg["time"], gg["PCR_Parasite_Prevalence"].rolling(14, min_periods=1).mean(), lw=1, label=f"seed {run}")
                ax.set_title(f"node {node}")
            axs[0, 0].legend(fontsize=7)
            fig.supxlabel("day")
            fig.supylabel("PCR prevalence (14-day mean)")
            fig.tight_layout()
            p = os.path.join(wdir, "seed_noise_prevalence.png")
            fig.savefig(p, dpi=120)
            print("Plot:", p)
        except Exception as e:
            print("(plot skipped:", e, ")")

    print("\nFiles in", wdir, ":", os.listdir(wdir))