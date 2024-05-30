# STARDUST
STARDUST: Spatio-Temporal Analysis of Regional Dynamics &amp; Unbiased Sorting of Transients

![STARDUST workflow](workflow.png)

STARDUST is a pipeline that captures intracellular calcium (Ca<sup>2+</sup>) dynamics in confined, local micro-domains across the territory of all astrocytes within the field of view. STARDUST is developed and maintained by the [Papouin lab](https://sites.wustl.edu/papouinlab/) at Washington Univeristy in St. Louis. We welcome suggestions and issue reporting via the Issue section on GitHub. 

This repository hosts code for the signal detection part of STARDUST. For a full description of the protocol, please refer to our paper [*STARDUST: a pipeline for the unbiased analysis of astrocyte regional calcium dynamics*](https://doi.org/10.1101/2024.04.04.588196). 

- [Overview](#overview)
- [Reference](#reference)
- [Updates](#updates)

# Overview
STARDUST builds upon [AQuA](https://github.com/yu-lab-vt/AQuA/), a popular open-source fluorescence analysis platform, to yield maps of regions of activity (ROAs) from patches of active signals, which can be combined with cell-segmentation and/or correlated to cellular morphology. Importantly, STARDUST makes no assumption regarding Ca<sup>2+</sup> propagation across ROAs, in line with the seemingly static nature of astrocyte Ca<sup>2+</sup> activity. Instead, STARDUST treats ROAs as independent units and focuses on decomposing Ca<sup>2+</sup> dynamics in a regionalized fashion, yielding around 20 ROAs per cell, or thousands across the field of view, and extracting fluorescence time-series, signals and signal features from each of them. A particular instantiation of the usefulness of STARDUST is in pharmacology experiments, where it can distinguish “stable ROAs” (active throughout the recording), from “ON ROAs” (inactive at baseline epoch and turned on during drug application), and “OFF ROAs” (active at baseline and turned off during drug application). Together, this makes STARDUST a user-friendly complement or alternative to the small number of publicly available algorithms and tools recently developed to tackle astrocyte Ca<sup>2+</sup> activity. This protocol includes step-by-step guidelines, tips on how to use STARDUST and outlines multiple output examples, providing a systematic walkthrough of its core functionalities, limitation, and tunability.

# Reference
Wu, Y., Dai, Y., Lefton, K. B., Holy, T. E. & Papouin, T. [*STARDUST: a pipeline for the unbiased analysis of astrocyte regional calcium dynamics.*](https://doi.org/10.1101/2024.04.04.588196) bioRxiv 2024.04.04.588196 (2024) doi:10.1101/2024.04.04.588196.

# Updates
2024/05/30: Make repository public. 