# Command line segmentation for multimodal/multichannel 3D volume data

Command line running a multimodal/channel image segmentation, based on:
* Principal Component Analysis (PCA), using itk library (https://www.itk.org)
* Support Vector Machines (SVM), using libsvm (https://www.csie.ntu.edu.tw/~cjlin/libsvm)
* Extended Random Walker (ERW), using vtk (http://www.vtk.org/) /eigen (http://eigen.tuxfamily.org)

## Building

svmerw requires CMake (https://cmake.org/) for setting up the build environment,
as well as the libraries mentioned above.
Download and build above libraries (see e.g. [the build instructions from our open_iA tool](https://github.com/3dct/open_iA/wiki/Building-open_iA)).

## Syntax
```
svmerw Mode (SeedFileName Flags MaskParams OutputFileName {InputFileName} AlgorithmParameters [GADParameters] | InputDirectory)
```
with 
```
Mode = "RW" | "ERW" (| "Norm")
Flag = "true" | "false"
```

The different `Mode` values stand for:
* `RW` for the Random Walker algorithm
* `ERW` for a segmentation pipeline employing PCA, SVM and ERW
* `Norm` for a normalization of the input

Modes `RW` and `ERW` share the following settings:
* `OutputFileName` the output file. Currently supported is only the MetaImage format (*.mha / *.mhd)
* `{InputFileName}` a list of input files; an arbitrary number of files can be specified, currently in either MetaImage (*.mha / *.mhd) or the custom Volume Stack (*.volstack) format. The first parameter that is interpretable as double is considered to be the first from the `AlgorithmParameters` set. Note that svmerw supports both both multimodal and multichannel data. Every input is considered a separate modality; if the input is given as Volume Stack, each volume in that stack is considered one channel of that modality. The total number of input channels is the sum of the number of channels for each modality (each modality can have a different number of channels).
* `SeedFile` the filename of an XML file with a format like in this example:
```
<?xml version="1.0" encoding="UTF-8"?>
<Labels>
    <Label id="0" name="Background">
        <Seed x="1" y="1" z="0"/>
	</Label>
	<Label id="1" name="Object1">
	    <Seed x="10" y="1" z="0"/>
	</Label>
	<Label id="2" name="Object2">
	    <Seed x="1" y="4" z="4"/>
	</Label>
</Labels>
```
Here, the pixel values from the input images will be considered.

Optionally, for the ERW pipeline, per Seed a value for each input channel can be specified, such as:
``` 
        <Seed x="1" y="1" z="0">
            <Value modality="0" component="0" value="2984"/>
            <Value modality="0" component="1" value="3880"/>
            <Value modality="1" component="0" value="4062"/>
		</Seed>
```
When seeds are specified like this, the ERW pipeline will not use the given x,y,z indices to look values in the input images, instead it will directly consider the values specified.

Note that the RW mode always requires seed coordinates in the given input image and will therefore not consider such given values.

### Mode `RW`

* `Flags = StoreRWProb UseGAD`
  Each of those is a Flag as defined above (taking a value of "true" or "false")
* `AlgorithmParameters = RWBeta {ModalityParameters}`
  * `RWBeta` determines the power in the Gaussian kernel used for distance calculation in the Random Walker, and typically takes values between 0.00001 and 10000
* `ModalityParameters = PCAComponents Weight DistanceFunction`
  These parameters need to be specified for each modality (that is, input file) given.
  * `PCAComponents` is the number of PCA components to consider. Valid values: 1 to the number of channels for that modality.
  * `Weight` the weight given to this particular modality
  * `DistanceFunction = "l1" | "l2" | "linf" | "sq" | "cos" | "js" | "kl" | "cs" | "em"
    determines the distance function used in the extended random walker for that modality:
	* l1: "Manhatten" / "Taxicab" Distance
	* l2: Euclidean Distance
	* linf: L-Infinity Norm (maximum distance of single component)
	* sq: the original distance function from Random walks (Grady [1]).
	* cos: cosine similarity (only meaningful for multichannel data)
	* js: Jensen-Shannon divergence (only meaningful for multichannel data)
	* kl: Kullback-Leibler divergence (only meaningful for multichannel data)
	* cs: Chi-Square distance (only meaningful for multichannel data)
	* em: Earth Mover's Distance

	
### Mode `ERW`

* `Flags = StoreSVMProb StoreRWProb UseGAD`
  Each of those is a Flag as defined above (taking a value of "true" or "false")
* `AlgorithmParameters =  ERWBeta ERWGamma ERWMaxIter SVMC SVMGamma SVMChannels ModalityParameters`
  * for `ERWBeta`, see description of `RWBeta` in mode `RW`
  * `ERWGamma` determines the weight of the prior model in comparison to the edge weights in the neighbourhood graph, useful values are in the range of 0.01 to 100.
  * `ERWMaxIter` limits the number of iterations that are performed during iterative solving of linear equations in the extended random walker, useful values range from 100 to 10000.
  * `SVMC` determines the penalty for misclassification. Smaller values try to maximize the separation, even if it means that more values are misclassified. Values found useful so far range from 0.01 to 10000.
  * `SVMGamma` is the free parameter for the radial basis function used as kernel in SVM. Values found useful so far range from 1e-12 to 10.
  * `SVMChannels` the number of input channels to consider in the SVM. Should be at least one, and below the number of  input channels over all modalities.
* `ModalityParameters` are defined the same as for `RW` mode


### Masking

`MaskParams = UseMask [MaskFileName]`
`UseMask = Flag`

You can specify a foreground mask (see UseMask / MaskFileName parameters above).
If UseMask is "true", then a MaskFileName parameter will be expected.
The volume given by that filename is expected to be binary and of the same dimension as the other input datasets,
and all parts of the output where the mask is 0 are also set to 0 in the output.


### Gradient Anisotropic Diffusion (GAD)

`GADParameters = Iterations Step Conductance`

As preprocessing filter to reduce noise, gradient anisotropic diffusion can be applied. 
For an explanation of the parameters please refer to the documentation of the used [GradientAnisotropicDiffusionImageFilter](https://itk.org/Doxygen/html/classitk_1_1GradientAnisotropicDiffusionImageFilter.html) by ITK.


## Examples:

To run the ERW pipeline on a simple, single-modality, single-channel dataset, with no probability output:

> svmerw ERW "seeds.xml" false false false false output.mhd input.mhd 10 1 1000 1 0.00001 1 1 1 sq


### Mode `Normalize`

The only parameter this mode requires is an input directory.
The (*.mhd) image files starting with "prob" in that directory will be "normalized", that is, in each pixel,
a factor is applied such that the sum of the images in that pixel is 1.


## References

[1] L. Grady. Random walks for image segmentation. IEEE Trans. Pattern Analysis and Machine Intelligence, 28(11):1768–1783, Nov. 2006. doi: [10.1109/TPAMI.2006.233](https://doi.org/10.1109/TPAMI.2006.233)