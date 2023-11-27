# `pSPSS`: Phenomenological symmetry protected seesaw scenario

Authors: Stefan Antusch, Jan Hajer, Bruno Oliveira, Johannes Rosskopp

The model file has been published at [feynrules.irmp.ucl.ac.be/wiki/pSPSS](https://feynrules.irmp.ucl.ac.be/wiki/pSPSS)

The motivation and implementation of the `pSPSS` is discussed in:
> Stefan Antusch, Jan Hajer, Johannes Rosskopp '**[Simulating lepton number violation induced by heavy neutrino-antineutrino oscillations at colliders](https://inspirehep.net/literature/2167345)**' e-Print: [2210.10738 [hep-ph]](https://arxiv.org/abs/2210.10738).

We kindly ask to cite  cite this publication when using the `pSPSS` model file.

An example analysis at the LHC including the necessary statistical framework to claim discovery of oscillations is presented in:
> Stefan Antusch, Jan Hajer, Johannes Rosskopp '**[Beyond lepton number violation at the HL-LHC: Resolving heavy neutrino-antineutrino oscillations](https://inspirehep.net/literature/2606505)**' e-Print: [2212.00562 [hep-ph]](https://arxiv.org/abs/2212.00562).

An initial study at the FCC-ee is presented in:
> Stefan Antusch, Jan Hajer, Bruno Oliveira '**[Heavy neutrino-antineutrino oscillations at the FCC-ee](https://inspirehep.net/literature/2687846)**' e-Print: [2308.07297 [hep-ph]](https://arxiv.org/abs/2308.07297).

A detailed calculation of the damping due to decoherence has been published in:
> Stefan Antusch, Jan Hajer, Johannes Rosskopp '**[Decoherence effects on lepton number violation from heavy neutrino-antineutrino oscillations](https://inspirehep.net/literature/2676341)** e-Print: [2307.06208 [hep-ph]](https://arxiv.org/abs/2307.06208).

## Model description 

The `pSPSS` describes the interactions of a pseudo-Dirac pair of two Majorana degrees of freedom $N_1$ and $N_2$ generically appearing in low-scale seesaw models protected by a lepton-number like symmetry.<br/>
In the lepton number conserving (LNC) limit the interactions of the symmetry protected seesaw scenario (SPSS) with the Standard Model are

$$\mathcal L_\text{SPSS}^L = \overline N_i i \cancel \partial N_i -  y_{\alpha1} \widetilde H^\dagger \bar \ell_\alpha N_1^c - \overline N_1 m_M^{} N_2^c + \text{H.c.}$$

The additional lepton number violating (LNV) interactions must be small in order to generate light neutrino masses.<br/>
Additionally, they introduce a small mass splitting $\Delta m$ between the two heavy neutrino mass eigenstates

$$m_{4/5}^{} = m_M^{} + \frac12 m_M^{} |\theta|^2 \mp \frac12 \Delta m$$

where the contribution with the active sterile mixing parameter $\theta = m_D / m_M$ with $m_D = y_1 v$ and $v\simeq174 \text{ GeV}$ appears already in the LNC case.<br/>
The smallness of the LNV interactions ensures unobservable collider effects, with the exception of heavy neutrino-antineutrino oscillations since these are a macroscopic interference effect.<br/>
At leading order, the oscillations between LNC and LNV processes induced by them can be described by 

$$P^{\text{LNC}/\text{LNV}}_\text{osc}(\tau) = (1 \pm \cos\left(\Delta m \tau \right) \exp(-\lambda))/2$$

where $\lambda$ captures potential damping of the oscillations due to decoherence.<br/>
Therefore, the details of the seesaw model besides the generated mass splitting can be neglected when simulating pseudo-Dirac heavy neutrinos as long as the neutrino-antineutrino oscillations are taken into account.

## `FeynRules` implementation 

The !FeynRules model file contains in addition to the Standard Model parameter as free parameter the heavy neutrino
* Majorana mass $m_M$ `Mmaj`
* mass splitting $\Delta M$ `deltaM`
* mixing parameter $\theta_\alpha$ `theta1`, `theta2`, `theta3`
* damping parameter $\lambda$ `damping`

### Antiparticles of neutral particles

The (anti-)particle character of the pseudo-Dirac heavy neutrinos is characterized by their interaction with (anti-)leptons.
A heavy neutrino is produced in association with an anti-lepton and decays into a lepton, while a heavy antineutrino is produced in association with a lepton and decays into an antilepton.
In order to extend this definition to interactions with light neutrinos it is necessary to define their (anti-)particle character as well.
On a conceptional level the definition for the light (anti-)neutrinos parallels the definition for the heavy (anti-)neutrinos.
In order to simulate interactions with definite light neutrino states we have implemented the light neutrinos as Dirac particles in the model file `pSPSS_Dirac_v1.0.fr`.

## `MadGraph` patch 

In order to generate events with heavy neutrino-antineutrino oscillations it is necessary to patch the `[pSPSS]/bin/internal/common_run_interface.py` file in !MadGraph.
For the LHC analysis we have replaced the original code 


```
for event in lhe:
        for particle in event:
            id = particle.pid
            width = param_card['decay'].get((abs(id),)).value
            if width:
                vtim = c * random.expovariate(width/cst)
                if vtim > threshold:
                    particle.vtim = vtim
        #write this modify event
        output.write(str(event))
    output.write('</LesHouchesEvents>\n')
    output.close()
```

with the modified code

```
mass_splitting = param_card.get_value('PSPSS', 2)
damping = param_card.get_value('PSPSS', 6)
for event in lhe:
    leptonnumber = 0
    write_event = True
    for particle in event:
        if particle.status == 1:
            if particle.pid in [11, 13, 15]:
                leptonnumber += 1
            elif particle.pid in [-11, -13, -15]:
                leptonnumber -= 1
    for particle in event:
        id = particle.pid
        width = param_card['decay'].get((abs(id),)).value
        if width:
            if id in [8000011, 8000012]:
                tau0 = random.expovariate(width / cst)
                if 0.5 * (1 + math.exp(-damping)*math.cos(mass_splitting * tau0 / cst)) >= random.random():
                    write_event = (leptonnumber == 0)
                else:
                    write_event = (leptonnumber != 0)
                vtim = tau0 * c
            else:
                vtim = c * random.expovariate(width / cst)
            if vtim > threshold:
                particle.vtim = vtim
    # write this modify event
    if write_event:
        output.write(str(event))
output.write('</LesHouchesEvents>\n')
output.close()
```

The complete code applicable to the Z-pole run of the FCC-ee reads

```
    def do_add_time_of_flight(self, line):
        print("Running patched do_add_time_of_flight")

        args = self.split_arg(line)
        #check the validity of the arguments and reformat args
        self.check_add_time_of_flight(args)

        event_path, threshold = args
        #gunzip the file
        if event_path.endswith('.gz'):
            need_zip = True
            misc.gunzip(event_path)
            event_path = event_path[:-3]
        else:
            need_zip = False

        import random
        try:
            import madgraph.various.lhe_parser as lhe_parser
        except:
            import internal.lhe_parser as lhe_parser

        logger.info('Add time of flight information on file %s' % event_path)
        lhe = lhe_parser.EventFile(event_path)
        output = open('%s_2vertex.lhe' % event_path, 'w')
        #write the banner to the output file
        output.write(lhe.banner)

        # get the associate param_card
        begin_param = lhe.banner.find('<slha>')
        end_param = lhe.banner.find('</slha>')
        param_card = lhe.banner[begin_param+6:end_param].split('\n')
        param_card = param_card_mod.ParamCard(param_card)

        cst = 6.58211915e-25 # hbar in GeV s
        c = 299792458000 # speed of light in mm/s
        sm_lepton_list = [11, 12, 13, 14, 15, 16]

        pspss_n_list = [8000011, 8000012]
        mass_splitting = param_card.get_value('PSPSS', 2)
        damping = param_card.get_value('PSPSS', 6)

        # Loop over all events
        for event in lhe:
            leptonnumber = 0
            for particle in event:
                if particle.status == 1:
                    if particle.pid in sm_lepton_list:
                        leptonnumber += 1
                    elif -particle.pid in sm_lepton_list:
                        leptonnumber -= 1

            write_event = True
            for particle in event:
                id = particle.pid
                width = param_card['decay'].get((abs(id),)).value

                if width:
                    if id in pspss_n_list:
                        tau0 = random.expovariate(width / cst)

                        if 0.5 * (1 + math.exp(-damping) * math.cos(mass_splitting * tau0 / cst)) >= random.random():
                            write_event = (leptonnumber == 0)
                        else:
                            write_event = (leptonnumber != 0)

                    else:
                        tau0 = random.expovariate(width / cst)

                    vtim = c * tau0
                    if vtim > threshold:
                        particle.vtim = vtim

            # write this modify event
            if write_event:
                output.write(str(event))

        output.write('</LesHouchesEvents>\n')
        output.close()
        files.mv('%s_2vertex.lhe' % event_path, event_path)

        if need_zip:
            misc.gzip(event_path)
```

## Workaround 

The small mass splitting can cause problems in the automatic calculation of the decay width.
One way to fix this problem is by replacing the argument of the square root of the return value of the function `calculate_apx_psarea` in the file `[MadGraph]/mg5decay/decay_objects.py` by its absolute value.
