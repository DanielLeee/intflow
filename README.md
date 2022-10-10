# Integrated Flow

This is the official implementation of the following paper:

> [**Integrated 3D Flow-based Multi-atlas Brain Structure Segmentation**](https://doi.org/10.1371/journal.pone.0270339)
> 
> Yeshu Li, Ziming Qiu, Xingyu Fan, Xianglong Liu, Eric I-Chao Chang, Yan Xu 
> 
> *[PLOS ONE](https://journals.plos.org/plosone/)*

This implementation is a generalized version of the initial one adopted in the paper. It is applicable to registration between two grayscale volumes in arbitrary dimensions (e.g., A\*B\*C vs D\*E\*F), tested on the complete MICCAI 2012 dataset.

## Quick Start

Specify any additional compilation flags in `Makefile` and run:

```shell
make
```

No libraries such as OpenCV are required.

Only the ANALYZE75 format is supported.

Registration:

```shell
./intFlow -fix fixImage -mov movImage -o flowFile energyFile \
          [-alpha alpha] \
          [-gamma gamma] \
          [-nlevels numLevels] \
          [-wsize windowSize] \
          [-niter numIterations] \
          [-nscale numScales] \
          [-grayavg averageGrayVal] \
          [-siftavg avarageSIFTVal]
```

Fusion:

```shell
./labelTransfer -train movImage1 movSeg1 movImage2 movSeg2 ... movImageK movSegK \
                -flow flowFile1 ... flowFileK \
                -test fixImage \
                -output predictedSegmentation \
                [-nlabels numLabels] \
                [-grayavg averageGrayVal] \
                [-siftavg avarageSIFTVal]
```

The default parameters are suggested for the original MICCAI 2012 dataset. For the ADNI dataset, a possible set of good parameters is
```shell
./intFlow -fix orig_AD001_L -mov orig_NL001_L -o trans_AD001_L_NL001_L.flow energy_AD001_L_NL001_L.txt -alpha 2 -nlevels 3 -wsize 2 -grayavg 70.0 -siftavg 1.0
```

An example of performing label fusion of three atlases:
```shell
./labelTransfer -train orig_NL001_L aseg_NL001_L orig_NL002_L aseg_NL002_L orig_NL003_L aseg_NL003_L -flow trans_AD001_L_NL001_L.flow trans_AD001_L_NL002_L.flow trans_AD001_L_NL003_L.flow -test orig_AD001_L -output latrseg_AD001_L -nlabels 2 -grayavg 70.0 -siftavg 1.0
```

Some data examples can be found [here](https://github.com/DanielLeee/intflow/releases/download/data/example_data.zip). The original code, obsolete though, can be found [here](https://github.com/DanielLeee/intflow/releases/download/data/replicate.zip) if you are interested in replicating the results deterministically.


## Citation

Please cite our work if you find it useful in your research:

```
@article{10.1371/journal.pone.0270339,
    doi = {10.1371/journal.pone.0270339},
    author = {Li, Yeshu AND Qiu, Ziming AND Fan, Xingyu AND Liu, Xianglong AND Chang, Eric I-Chao AND Xu, Yan},
    journal = {PLOS ONE},
    publisher = {Public Library of Science},
    title = {Integrated 3d flow-based multi-atlas brain structure segmentation},
    year = {2022},
    month = {08},
    volume = {17},
    url = {https://doi.org/10.1371/journal.pone.0270339},
    pages = {1-37},
    abstract = {MRI brain structure segmentation plays an important role in neuroimaging studies. Existing methods either spend much CPU time, require considerable annotated data, or fail in segmenting volumes with large deformation. In this paper, we develop a novel multi-atlas-based algorithm for 3D MRI brain structure segmentation. It consists of three modules: registration, atlas selection and label fusion. Both registration and label fusion leverage an integrated flow based on grayscale and SIFT features. We introduce an effective and efficient strategy for atlas selection by employing the accompanying energy generated in the registration step. A 3D sequential belief propagation method and a 3D coarse-to-fine flow matching approach are developed in both registration and label fusion modules. The proposed method is evaluated on five public datasets. The results show that it has the best performance in almost all the settings compared to competitive methods such as ANTs, Elastix, Learning to Rank and Joint Label Fusion. Moreover, our registration method is more than 7 times as efficient as that of ANTs SyN, while our label transfer method is 18 times faster than Joint Label Fusion in CPU time. The results on the ADNI dataset demonstrate that our method is applicable to image pairs that require a significant transformation in registration. The performance on a composite dataset suggests that our method succeeds in a cross-modality manner. The results of this study show that the integrated 3D flow-based method is effective and efficient for brain structure segmentation. It also demonstrates the power of SIFT features, multi-atlas segmentation and classical machine learning algorithms for a medical image analysis task. The experimental results on public datasets show the proposed method’s potential for general applicability in various brain structures and settings.},
    number = {8},
}
```

## Acknowledgement

A large part of the code is borrowed from [Label Transfer](http://people.csail.mit.edu/celiu/LabelTransfer).
