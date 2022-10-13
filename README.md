# COMPAS CRA

![build](https://github.com/blockresearchgroup/compas_cra/workflows/build/badge.svg)
[![GitHub - License](https://img.shields.io/github/license/blockresearchgroup/compas_cra.svg)](./LICENSE)
[![pip downloads](https://img.shields.io/pypi/dm/compas_cra)](https://pypi.python.org/project/compas_cra)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/compas_cra.svg)](https://pypi.python.org/project/compas_cra)
[![PyPI - Latest Release](https://img.shields.io/pypi/v/compas_cra.svg)](https://pypi.python.org/project/compas_cra)
[![DOI](https://zenodo.org/badge/374677757.svg)](https://zenodo.org/badge/latestdoi/374677757)

[Coupled Rigid-Block Analysis (CRA)](https://doi.org/10.1016/j.cad.2022.103216) implementation using [COMPAS](https://compas.dev/) framework.

> developed with <span style="color: #e25555;">&#9829;</span> by [Gene Ting-Chun Kao](https://geneatcg.com) 

To find out more about CRA, please refer to our paper in the CAD Computer-Aided Design journal: 
[https://doi.org/10.1016/j.cad.2022.103216](https://doi.org/10.1016/j.cad.2022.103216 ) 

#### Coupled Rigid-Block Analysis: Stability-Aware Design of Complex Discrete-Element Assemblies


![image](https://github.com/BlockResearchGroup/compas_cra/blob/main/docs/_images/cra_bridge.png?raw=true)

##### Abstract

The rigid-block equilibrium (RBE) method uses a penalty formulation to
measure structural infeasibility or to guide the design of stable
discrete-element assemblies from unstable geometry.
However, RBE is a purely force-based formulation,
and it incorrectly describes stability when
complex interface geometries are involved.
To overcome this issue, this paper introduces
the coupled rigid-block analysis (CRA) method,
a more robust approach building upon RBE’s strengths.
The CRA method combines equilibrium and kinematics in a penalty formulation
in a nonlinear programming problem.
An extensive benchmark campaign is used to show how CRA enables
accurate modelling of complex three-dimensional discrete-element assemblies
formed by rigid blocks.
In addition, an interactive stability-aware design process to
guide user design towards structurally-sound assemblies is proposed.
Finally, the potential of our method for real-world problems are demonstrated
by designing complex and scaffolding-free physical models.


##### Please cite our work if you use CRA in your research

###### Paper

```latex
@article{kao2022coupled,
    title     = {Coupled Rigid-Block Analysis: Stability-Aware Design of Complex Discrete-Element Assemblies},
    author    = {Kao, Gene Ting-Chun and Iannuzzo, Antonino and Thomaszewski, Bernhard and Coros, Stelian and Van Mele, Tom and Block, Philippe},
    journal   = {Computer-Aided Design},
    pages     = {103216},
    year      = {2022},
    publisher = {Elsevier},
    doi       = {10.1016/j.cad.2022.103216},
    url       = {https://doi.org/10.1016/j.cad.2022.103216}
}
```

###### Software implementation

```latex
@misc{compas-cra,
    title  = {{COMPAS CRA}: Coupled Rigid-Block Analysis ({CRA}) for the {COMPAS} framework},
    author = {Gene Ting-Chun Kao},
    note   = {https://github.com/BlockResearchGroup/compas{\_}cra},
    year   = {2020-2022},
    doi    = {10.5281/zenodo.7043135},
    url    = {https://doi.org/10.5281/zenodo.7043135},
}
```

##### Read the docs
[https://blockresearchgroup.github.io/compas_cra](https://blockresearchgroup.github.io/compas_cra)

##### Build the docs locally

   ```bash
   $ pip install -r requirements-dev.txt
   $ invoke docs
   $ open dist/docs/index.html  # or open index.html in compas_cra/dist/docs/
   ```

##### Examples to reproduce our paper results

See examples in [docs](https://blockresearchgroup.github.io/compas_cra/latest/examples.html) or try them in [docs/examples](https://github.com/BlockResearchGroup/compas_cra/blob/main/docs/examples). 
