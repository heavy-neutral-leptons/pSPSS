# `pSPSS`: Phenomenological symmetry protected seesaw scenario

The motivation and implementation of the `pSPSS` is discussed in:

Stefan Antusch, Jan Hajer, Johannes Rosskopp '**Simulating lepton number violation induced by heavy neutrino-antineutrino oscillations at colliders**' e-Print: [2210.10738 [hep-ph]](https://arxiv.org/abs/2210.10738).


Please cite this publication when you use the `pSPSS` model file.

## Model description

The `pSPSS` describes the interactions of a pseudo-Dirac pair of two Majorana degrees of freedom $N_1$ and $N_2$ generically appearing in low-scale seesaw models protected by a lepton-number like symmetry.<br/>
In the lepton number conserving (LNC) limit the interactions of the symmetry protected seesaw scenario (SPSS) with the Standard Model are

$$\mathcal L_\text{SPSS}^L = \overline N_i \, i \slashed \partial N_i -  y_{\alpha1} \widetilde H^\dagger \bar \ell_\alpha N_1^c - \overline N_1 m_M^{} N_2^c + \text{H.c.}$$

The additional lepton number violating (LNV) interactions must be small in order to generate light neutrino masses.<br/>
Additionally, they introduce a small mass splitting $\Delta m$ between the two heavy neutrino mass eigenstates

$$m_{4/5}^{} = m_M^{} + \frac12 m_M^{} |\theta|^2 \mp \frac12 \Delta m$$

where the contribution with the active sterile mixing parameter $\theta = m_D / m_M$ with $m_D = y_1 v$ and $v\simeq174 \text{ GeV}$ appears already in the LNC case.<br/>
The smallness of the LNV interactions ensures unobservable collider effects, with the exception of heavy neutrino-antineutrino oscillations since these are a macroscopic interference effect.<br/>
At leading order, the oscillations between LNC and LNV processes induced by them can be described by

$$P^{\text{LNC}/\text{LNV}}_\text{osc}(\tau) = (1 \pm \cos\left(\Delta m \tau \right) \exp(-\lambda))/2$$

where $\lambda$ captures potential damping of the oscillations due to decoherence.<br/>
Therefore, the details of the seesaw model besides the generated mass splitting can be neglected when simulating pseudo-Dirac heavy neutrinos as long as the neutrino-antineutrino oscillations are taken into account.

## FeynRules implementation

The FeynRules model file contains in addition to the Standard Model parameter as free parameter the heavy neutrino

* Majorana mass $m_M$ `Mmaj`
* mass splitting $\Delta M$ `deltaM`
* mixing parameter $\theta_\alpha$ `theta1`, `theta2`, `theta3`
* damping parameter $\lambda$ `damping`

## MadGraph patch

In order to generate events with heavy neutrino-antineutrino oscillations it is necessary to patch the `[pSPSS]/bin/internal/common_run_interface.py` file in !MadGraph using


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
                if 0.5 * (1 + math.exp(=damping)*math.cos(mass_splitting * tau0 / cst)) >= random.random():
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

which must replace the original code

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

The small mass splitting can cause problems in the automatic calculation of the decay width.
One way to fix this problem is by replacing the argument of the square root of the return value of the function `calculate_apx_psarea` in the file `[MadGraph]/mg5decay/decay_objects.py` by its absolute value.
