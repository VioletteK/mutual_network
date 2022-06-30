{

    'run_time': 3000, # ms
    'dt': 0.1, # ms

    'Populations' : {
        'drive' : {
            #'n' : 25*25,
            'n' : 100*100,
            'type': sim.SpikeSourcePoisson,
            'cellparams' : {
                'start':0.0,
                'rate':4.,
                'duration': 3000.0
            }
        },

       'py' : {
            'n': 100*100, # units
            'type': sim.EIF_cond_alpha_isfa_ista,
            'structure' : Grid2D(aspect_ratio=1, dx=1.0, dy=1.0, fill_order='sequential', rng=sim.NumpyRNG(seed=2**32-1)),
            'cellparams': {
                'e_rev_I'    : -80,   # mV, reversal potential of excitatory synapses
                'e_rev_E'    : 0,     # mV, reversal potential of inhibitory synapses
                'tau_syn_E'  : 5.0,   # ms, time constant of excitatory synaptic short-term plasticity, YgerBoustaniDestexheFregnac2011
                'tau_syn_I'  : 5.0,   # ms, time constant of excitatory synaptic short-term plasticity, YgerBoustaniDestexheFregnac2011
                'tau_refrac' : 5.0,   # ms, refractory period
                'v_reset'    : -65.0, # mV, reset after spike
                'v_thresh'   : -50.0, # mV, spike threshold (modified by adaptation)
                'delta_T'    : 2.,    # mV, steepness of exponential approach to threshold
                'cm'         : 0.150, # nF, tot membrane capacitance
                'a'          : 4.,    # nS, conductance of adaptation variable
                'tau_m'      : 15.0,  # ms, time constant of leak conductance (cm/gl)
                'v_rest'     : -65.0, # mV, resting potential E_leak
                'tau_w'      : 500.0, # ms, time constant of adaptation variable
                'b'          : .02,   # nA, increment to adaptation variable
            },
        },
        'inh' : {
            'n': 50*50, #{'ref':'py','ratio':0.25},
            'type': sim.EIF_cond_alpha_isfa_ista,
            'structure' : Grid2D(aspect_ratio=1, dx=2.0, dy=2.0, fill_order='sequential', rng=sim.NumpyRNG(seed=2**32-1)),
            'cellparams': {
                'e_rev_I'    : -80,   # mV, reversal potential of excitatory synapses
                'e_rev_E'    : 0,     # mV, reversal potential of inhibitory synapses
                'tau_syn_E'  : 5.0,   # ms, time constant of excitatory synaptic short-term plasticity, YgerBoustaniDestexheFregnac2011
                'tau_syn_I'  : 5.0,   # ms, time constant of inhibitory synaptic short-term plasticity, YgerBoustaniDestexheFregnac2011
                'tau_refrac' : 5.0,   # ms, refractory period
                'v_reset'    : -65.0, # mV, reset after spike
                'v_thresh'   : -50.0, # mV, spike threshold (modified by adaptation)
                'delta_T'    : 0.5,   # mV, steepness of exponential approach to threshold
                'cm'         : 0.150, # nF, tot membrane capacitance
                'a'          : 0.0,   # nS, conductance of adaptation variable
                'tau_m'      : 15.0,  # ms, time constant of leak conductance (cm/gl)
                'v_rest'     : -65.0, # mV, resting potential E_leak
                'tau_w'      : 500.0, # ms, time constant of adaptation variable
                'b'          : 0.0,   # nA, increment to adaptation variable
            },
        },
    },

    'Projections' : {
        'drive_py' : {
            'source' : 'drive',
            'target' : 'py',
            'space' :  sim.Space(periodic_boundaries=((0,100), (0,100), None)), # torus
            'connector' : sim.FixedProbabilityConnector(.01, rng=sim.random.NumpyRNG(2**32-1)),
            'synapse_type' : sim.StaticSynapse(),
            # 'weight' : .003, # uS # 25*25 *1000 *.008 = 5000
            'weight' : .0005, # uS # 100*100 *1000 *0.005 = 5000
            'receptor_type' : 'excitatory',
            'save_connections':False,
            'print_statistics':False,
        },
        'drive_inh' : {
            'source' : 'drive',
            'target' : 'inh',
            'space' :  sim.Space(periodic_boundaries=((0,100), (0,100), None)), # torus
            'connector' : sim.FixedProbabilityConnector(.01, rng=sim.random.NumpyRNG(2**32-1)),
            'synapse_type' : sim.StaticSynapse(),
            # 'weight' : .003, # uS # 25*25 *1000 *.008 = 5000
            'weight' : {'ref':'drive_py'}, # uS # 100*100 *1000 *0.005 = 5000
            'receptor_type' : 'excitatory',
            'save_connections':False,
            'print_statistics':False,
        },

        'py_py' : {
            'source' : 'py',
            'target' : 'py',
            'space' :  sim.Space(periodic_boundaries=((0,100), (0,100), None)), # torus
            'connector' : sim.DistanceDependentProbabilityConnector("14*exp(-1.2*d)", allow_self_connections=False, rng=sim.NumpyRNG(2**32-1)), # radius 300um
            'weight' : .001, # uS
            'synapse_type' : sim.StaticSynapse(),
            'delay' : .5, # ms
            'receptor_type' : 'excitatory',
            'save_connections':False,
            'print_statistics':False,
        },
        'py_inh' : {
            'source' : 'py',
            'target' : 'inh',
            'space' :  sim.Space(periodic_boundaries=((0,100), (0,100), None)), # torus
            'connector' : sim.DistanceDependentProbabilityConnector("24*exp(-1.5*d)", allow_self_connections=False, rng=sim.NumpyRNG(2**32-1)), # radius 100um
            'weight' : .001, # uS
            'synapse_type' : sim.StaticSynapse(),
            'delay' : .5, # ms,
            'receptor_type' : 'excitatory',
            'save_connections':False,
            'print_statistics':False,
        },
        'inh_inh' : {
            'source' : 'inh',
            'target' : 'inh',
            'space' :  sim.Space(periodic_boundaries=((0,100), (0,100), None)), # torus
            'connector' : sim.DistanceDependentProbabilityConnector("14*exp(-1.2*d)", allow_self_connections=False, rng=sim.NumpyRNG(2**32-1)), # radius 300um
            'weight' : .005, # uS
            'synapse_type' : sim.StaticSynapse(),
            'delay' : .5, # ms,
            'receptor_type' : 'inhibitory',
            'save_connections':False,
            'print_statistics':False,
        },
        'inh_py' : {
            'source' : 'inh',
            'target' : 'py',
            'space' :  sim.Space(periodic_boundaries=((0,100), (0,100), None)), # torus
            'connector' : sim.DistanceDependentProbabilityConnector("24*exp(-1.5*d)", allow_self_connections=False, rng=sim.NumpyRNG(2**32-1)), # radius 200um
            'weight' : .005, # uS
            'synapse_type' : sim.StaticSynapse(),
            'delay' : .5, # ms,
            'receptor_type' : 'inhibitory',
            'save_connections':False,
            'print_statistics':False,
        },

    },


    'Recorders' : {
        'py' : {
            'spikes' : 'all',
            'v' : {
                'MUA': True,
                'x': 18, #Correspond à la limite inférieure gauche du carré de 10x10 centré en 32 donc 32-5=27
                'y': 18,
                'size': 64,
            },
            'gsyn_exc' : {
                'start' : 0,
                'end' : 10,
            },
            'gsyn_inh' : {
                'start' : 0,
                'end' : 10,
            },
        },
        'inh' : {
            'spikes' : 'all',
            # 'v' : {
            #     'start' : 0,
            #     'end' : 100,
            # },
            'gsyn_exc' : {
                'start' : 0,
                'end' : 10,
            },
            'gsyn_inh' : {
                'start' : 0,
                'end' : 10,
            },
        },
    },


    'Modifiers' :{
    },

    'Injections' : {
        'py' : { # 'modKey':{modVal}
             'source' : sim.StepCurrentSource,
             'amplitude' : [.4, .0], # default
             'start' : [1500., 2000.], # long duration
             'stop' : 0.0,
             #'cellidx' : [27,28,35,36],
             #'cellidx' : 50+(100*50), # On prend la cellule au milieu du carré donc on monte 50 lignes à partir d'en bas (sur un tableau de 100 cellules) puis on se décale à la 50eme colonne
             #'cellidx' : 32+(64*32), # On prend la cellule au milieu du carré donc on monte 32 lignes à partir d'en bas (sur un tableau de 64 cellules) puis on se décale à la 32eme colonne
             #'cellidx' : [2015,2016,2017,2079,2080,2081,2143,2144,2145], #3x3 cellinjected
             #[2015, 48], [2016, 49], [2017, 50],  [2079, 59], [2080, 60], [2081, 61], , [2143, 70], [2144, 71], [2145, 72],
             #'cellidx' : [5149,5150,5151,5049,5050,5051,4949,4950,4951], #3x3 cellinjected

             ##### 7x7 ceEll injectedD
             #'cellidx' : [5347,5348,5349,5350,5351,5352,5353,5247,5248,5249,5250,5251,5252,5253,5147,5148,5149,5150,5151,5152,5153,5047,5048,5049,5050,5051,5052,5053,4947,4948,4949,4950,4951,4952,4953,4847,4848,4849,4850,4851,4852,4853,4747,4748,4749,4750,4751,4752,4753]
             ##### 5x5 cell injected
             # 'cellidx' : [5248,5249,5250,5251,5252,5148,5149,5150,5151,5152,5048,5049,5050,5051,5052,4948,4949,4950,4951,4952,4848,4849,4850,4851,4852],
             #


             # 'cellidx' : 32+(64*32), # On prend la cellule au milieu du carré donc on monte 32 lignes à partir d'en bas (sur un tableau de 64 cellules) puis on se décale à la 32eme colonne
             # 'cellidx' : [32*32, 7543, 6536], # On prend la cellule au milieu du carré
             ###### 20x20 cell injected :
             'cellidx' : [i*100 + j for i in range(30,70) for j in range(30,70)]
         },
    },

    'Analysis' : {
        # 'subsampling': 1000, # number of randomly selected units spiketrains for analysis (10% of the total as for gentic or rabies labelling)

        'scores' : ['py'],
        # 'scores' : ['py','inh'],

        'transient' : 2, # ms

        'Vm' : {
            'py'
        },

         'ConductanceBalance' : {
            'py':{
                'trials': ['default'], # for which trials the analysis has to be computed
            },
            # 'inh':{
            #     'trials': ['default'], # for which trials the analysis has to be computed
            # },
        },

        'FiringRate' : {
            'bin': 10, # ms
            'py':{
                'firing': [0,50],
            },
            'inh':{
                'firing': [0,50],
            },
        },
        #  'Rasterplot' : {
        #      'py':{
        #         'limits': [(0,63),(0,63)], # coords: [(from x, to x), (from y, to y)]
        #          'color': 'red',
        #       },
        #       'inh':{
        #          'limits': [(0,63),(0,63)], # coords: [(from x, to x), (from y, to y)]
        #          'color': 'blue',
        #     },
        #     'type': '.png',
        # 'type': '.svg',
        #     'interval': False, # all
        #     # 'interval': [2000.,3000.], # ms # from 2s to 3s
        #     'dpi':800,
        # },
    },

 }
