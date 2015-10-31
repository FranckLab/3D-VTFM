This repositry contains the Matlab m-files to run our LD 3D TFM algorithm. The first part of the package includes the FIDVC algorithm, which is how the 3D displacement fields are calculated. The second part of the code package includes converting those displacement fields into 3D surface tractions as described in [Toyjanova, Hannen et al., Soft Matter, 2014](http://pubs.rsc.org/en/content/articlelanding/2014/sm/c4sm01271b#!divAbstract).

Most of the Matlab file structures are very similar to our LD 3D TFM algorithm. The current code structure is set up to compute viscoelastic tractions based on the constitutive model depicted in our Soft Matter paper ([Toyjanova, Hannen et al., Soft Matter, 2014](http://pubs.rsc.org/en/content/articlelanding/2014/sm/c4sm01271b#!divAbstract)), but the user can simply tune the code to fit the specifics of their mechanical viscoelastic model. Feel free to contact us with questions or with help on implementing a user-specific model into our 3D VTFM code structure. 

### Important pages
* [Download latest version v1.0!](https://github.com/FranckLab/3D-VTFM/releases)
* [Example data](https://drive.google.com/folderview?id=0ByhZqlrbo5srQ1RDQ3pvZjVRSUU&usp=sharing)
* [FAQ](https://github.com/FranckLab/3D-VTFM#faq)
* [Questions/Issues](https://github.com/FranckLab/3D-VTFM/issues)
* [Bug Fixes/history](https://github.com/FranckLab/3D-VTFM/wiki/Bug-Fixes!)
* [Cite](https://github.com/FranckLab/3D-VTFM#cite)
* [Franck Lab](http://franck.engin.brown.edu)
 
## Running 3D-VTFM

### C Compiler
To run you need a compatible C compiler. Please see
(http://www.mathworks.com/support/compilers/R2015a/index.html)

### Input 3D Image Stack Requirements
* To check if the 3D image stack have the required speckle pattern and intensity values for correlation please use our [DVC simulator](https://github.com/FranckLab/DVC-Simulator).
* The 3D image stack need to be saved in a 3 dimensional matrix (intensity values are stored at x, y and z position) in **vol*.mat** files.  
* We recommend that the input image stack at each dimension should have at least 1.5 times of the subset size as the number of pixels. The default subset size is 128x128x64, so we recommend that the minimum input volume size should be 192x192x96.
* The size of the input image stack should be divisible by 0.5 times the size of the subset. 

### Running included example case
1. Make sure that the main files and the supplemental m files (from file exchange) are added to the path in Matlab.
2. Download and save the [example volume data](https://drive.google.com/folderview?id=0ByhZqlrbo5srQ1RDQ3pvZjVRSUU&usp=sharing) in the example folder. 
3. Run the exampleRunFile.m file

### Health warning!
FIDVC in the 3D-VTFM requires a 3D stack to be read in, which depending on the volume size can require a **large amount of RAM** in Matlab.

## Files
* Main files
 - calculateNormals.m
 - calculateSurfaceUi.m
 - findSurface.m
 - fun3VDTFM.m
 - removeOutliers.m

* Supplement m files from the MATLAB file exchange:
 - gridfit.m
 - inpaint_nans.m
 - inpaint_nans3.m

* Example Run files
 - exampleRunFile.m
 - [example volume data](https://drive.google.com/folderview?id=0ByhZqlrbo5srQ1RDQ3pvZjVRSUU&usp=sharing) (vol00.mat, vol01.mat).

## FAQ
**What are the requirements for the input 3D image stack?**

Please refer to [input 3D Image Stack Requirements](https://github.com/FranckLab/FIDVC/blob/master/README.md#input-3d-image-stack-requirements).


## Cite
If used please cite:
[Toyjanova J., Hannen, E., Bar-Kochba E., Darling, E.M., Henann, D.L., and Franck, C., (2014) 3D Viscoelastic Traction Force Microscopy. Soft Matter doi: 10.1039/c4sm01271b](http://pubs.rsc.org/en/content/articlelanding/2014/sm/c4sm01271b#!divAbstract)

```bibtex
@article{toyjanova20143d,
  title={3D Viscoelastic traction force microscopy},
  author={Toyjanova, Jennet and Hannen, Erin and Bar-Kochba, Eyal and Darling, Eric M and Henann, David L and Franck, Christian},
  journal={Soft matter},
  volume={10},
  number={40},
  pages={8095--8106},
  year={2014},
  publisher={Royal Society of Chemistry}
}
```

## Contact and support
For questions, please first refer to [FAQ](https://github.com/FranckLab/3D-VTFM#faq) and [Questions/Issues](https://github.com/FranckLab/3D-VTFM/issues). Add a new question if similar issue hasn't been reported. We shall help you at the earliest. The author's contact information can be found at [Franck Lab](http://franck.engin.brown.edu).
