Using PlottingTools/PlotWithVarial
----------------------------------

In PlottingTools/PlotWithVarial you can find a plotting tool that runs in parallel and is configured in python. It creates a webpage on which the plots are shown and the plots in different categories are linked. These tools rely on Varial, see https://github.com/HeinerTholen/Varial. Its installation and updating is automatically done when you ``source env.sh`` in the AnalysisTools base directory.

The plotting tool is configured with a config file, where a couple of variables are defined. This is an example:

    name = 'VHbbPlots'
    input_pattern = '/some/path/V25_VHbb_runonskim_20180212/haddjobs/sum_%s.root'
    weight = 'weight*puWeight*sign(genWeight)'
    enable_reuse_step = True  # try to find output on disk and don't run a step if present
    
    from main_samples import the_samples_dict, sample_colors
    from main_selections import the_category_dict

The last two lines import the samples and category definitions. Have a look at these files to see how ``the_category_dict`` and ``the_samples_dict`` are defined. This exact structure is expected by the tool.

IMPORTANT: ``input_pattern`` must contain ``%s``! This is where the input token from the sample definition is inserted. Wildcards are allowed.

Once Varial is installed, alter ``input_pattern`` in the config file an run this script. For the first time running it might be good to start with ``demo_config.py``, which is identical to ``main_config.py``, but only runs two samples in the signal region:

    ./run_plot_from_tree.py demo_config.py

