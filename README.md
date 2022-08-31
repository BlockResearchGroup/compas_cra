# COMPAS CRA

[Coupled Rigid-Block Analysis (CRA)](https://doi.org/10.1016/j.cad.2022.103216) implementation using [COMPAS](https://compas.dev/) framework.

> developed with <span style="color: #e25555;">&#9829;</span> by [Gene Ting-Chun Kao](https://geneatcg.com) 

To find out more about CRA, please refer to our paper in the CAD Computer-Aided Design journal: 
[https://doi.org/10.1016/j.cad.2022.103216](https://doi.org/10.1016/j.cad.2022.103216 ) 

#### Coupled Rigid-Block Analysis: Stability-Aware Design of Complex Discrete-Element Assemblies


![image](./docs/_images/cra_bridge.png)

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

```latex
@article{kao2022coupled,
  title={Coupled Rigid-Block Analysis: Stability-Aware Design of Complex Discrete-Element Assemblies},
  author={Kao, Gene Ting-Chun and Iannuzzo, Antonino and Thomaszewski, Bernhard and Coros, Stelian and Van Mele, Tom and Block, Philippe},
  journal={Computer-Aided Design},
  pages={103216},
  year={2022},
  publisher={Elsevier}
}
```

Read the docs: [https://blockresearchgroup.github.io/compas_cra](https://blockresearchgroup.github.io/compas_cra)

Build the docs locally: 

   ```bash
   $ pip install -r requirements-dev.txt
   $ invoke docs
   $ open dist/docs/index.html  # or open index.html in compas_cra/dist/docs/
   ```

