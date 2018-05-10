# Pose Changes From a Different Point of View
This article introduces the Pose Change Group(PCG) for viewing changes in pose from a space-fixed reference frame. In group-theoretic terms, the space of pose changes can be endowed with a direct product structure, with bi-invariant metric. The work has been published in the ASME Journal of Mechanism and Robotics [[pdf]](https://rpk.lcsr.jhu.edu/wp-content/uploads/2018/05/kinematics-new-view-for-publish-final.pdf).

Authors: Gregory S. Chirikjian (<gchirik1@jhu.edu>), Robert Mahony (<robert.mahony@anu.edu.au>), Sipu Ruan (<ruansp@jhu.edu>), Jochen Trumpf (<jochen.trumpf@anu.edu.au>)

## Introduction
For more than a century, rigid-body displacements have been viewed as affine transformations described as homogeneous transformation matrices wherein the linear part is a rotation matrix. In group-theoretic terms, this classical description makes rigid-body motions a semi-direct product. The distinction between a rigid-body displacement of Euclidean space and a change in pose from one reference frame to another is usually not articulated well in the literature. Here we show that, remarkably, when changes in pose are viewed from a space-fixed reference frame, the space of pose changes can be endowed with a direct product group structure, which is different than the semi-direct product structure of the space of motions. We then show how this new perspective can be applied more naturally to problems such as monitoring the state of aerial vehicles from the ground, or the cameras in a humanoid robot observing pose changes of its hands.

## Code Description
This repository publishes the source code that generates the applications for the idea of a Pose Change Group(PCG), with comparisons with Special Euclidean Group(SE):
1. Path generations in 2D and 3D cases, testing files are "pathGen2D.m" and "pathGen3D.m";
2. Trajectory interpolations in 3D, testing file is "path_multiFrames_3D.m".

The key functions for interpolation between two end points and multiple middle poins are "interpX.m" and "interpMultiPt.m" respectively, both of which implement the interpolations in Euclidean space(R), and Lie groups(SE and PCG).

## Citations
To reference the Pose Change Group in publications, please cite:

<cite>Chirikjian, G.S., Mahony, R., Ruan, S. and Trumpf, J., 2018. Pose Changes From a Different Point of View. Journal of Mechanisms and Robotics, 10(2), p.021008.</cite>

```
@article{chirikjian2018pose,
  title={Pose Changes From a Different Point of View},
  author={Chirikjian, Gregory S and Mahony, Robert and Ruan, Sipu and Trumpf, Jochen},
  journal={Journal of Mechanisms and Robotics},
  volume={10},
  number={2},
  pages={021008},
  year={2018},
  publisher={American Society of Mechanical Engineers}
}
```
