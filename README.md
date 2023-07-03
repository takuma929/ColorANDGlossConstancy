# Repository for color and gloss constancy project

A repository to store raw experimental data, images, and MATLAB source codes to reproduce result figures in an associated publication.


Takuma Morimoto, Arash Akbarinia, Katherine Storrs, Jacob R. Cheeseman, Hannah E. Smithson, Karl R. Gegenfurtner and Roland W. Fleming, “Color and gloss constancy under diverse lighting environments”. bioRxiv. https://doi.org/10.1101/2022.12.09.519756 (in press, Journal of Vision).

![alt text](https://github.com/takuma929/ColorANDGlossConstancy/blob/main/thumbnail.png?raw=true)


**How to Use** 

Change directory to the main repository first, and then run "main.m" which will generate all images. All images will be save to figs folder.

Make sure that "Psychtoolbox 3" is installed.

**Files**

There are 4 main codes to generate figures:

"main.m" to run other three codes

"CGC_figParameters.m" to define parameters for figure generation (e.g. fontsize, fontname)

"CGC_fig_Exp1.m" to generate figures in expeirment 1

"CGC_fig_Exp2.m" to generate figures in expeirment 2

The "data" folder stores image data (imgs folder) and raw experimental data (rawdata folder).
There are also imgStats_exp1.mat and imgStats_exp2.mat which store pre-computed image statistics for images used in experiments 1 and 2.

**Raw data**

"exp 1" folder stores each individiaul data for Experiment 1.
The filename format is "{subjectName}_session{session number}.mat".
Each data stores 2 variables.
(1) "responses" stores human judgement data for hue, lightness, chroma and Pellacini_c for 38 images (36 test images + 2 control condition).
(2) "groundtruth" stores associated ground-truth value for each image, and imageN indicates test images (1-36) and control images (control image 1 = -2, control image 2 = -1).

Image numbers 1 - 36 corresponds to test images in the order of "img1_natural.mat", "img1_gamutrotated.mat", "img1_phasescrambled.mat", "img2_natural.mat", ... "img12_phasescrambled.mat". 

"exp 2" folder stores each individiaul data for Experiment 2.
The filename format is "{subjectName}_session{session number}.mat".
Each data stores 2 variables: 
(1) "responses" stores human judgement data for Pellacini_c for 218 images (216 test images + 2 control condition). 
(2) "groundtruth" stores associated ground-truth value for each image, and imageN indicates test images (1-216) and control images (control image 1 = -2, control image 2 = -1).

Image numbers 1 - 216 corresponds to test images in the order of "img_lightprobe1_obj1.mat", "img_lightprobe1_obj2.mat", ..., "img_lightprobe12_obj1.mat",...,"img_lightprobe12_obj18.mat".


**Image data**

"exp1" stores 36 test images. Each mat file stores XYZ1931, srgb, and alpha images.
"exp2" stores 216 test images. Each mat file stores XYZ1931, and srgb images.

"lightprobe" folder stores thumbnail images for 12 lightprobes we used in the experiment.
