Feature Selection Multiple Kernel Learning

This software was developed by
Jose A. Seoane, Ian Day, Tom Gaunt and Colin Campbell

This software is based in SimpleMKL, developed by A. Rakotomamonjy, F. Bach, Y. Grandvalet, S. Canu 
(SimpleMKL,  Journal of Machine Learning Research, Vol. 9, pp 2491-2521, 2008)


The code was originally developed to be used with METABRIC dataset.

The data created for this example (randomData.mat) was created artificially (generateData.m) simulating the original distribution of the metabric dataset. 

The pathways used in this example are KEGG pathways and are stored in pathways.mat

This code runs in Matlab 7.12 with Bioinformatics Toolbox

The entry file is: PathwaySelectionFSMKL_main.m

INSTRUCTIONS:
The PathwaySelectionFSMKL.m example run the MKL two times for survival prediction in clinical, EXP and CNV datasets. First time with all data, and the second time without the outliers that selected in the first time based in its probability.

input of the algorithm:
- Input Datasets: The algorithm can work with one or several datasets. We have probed that most of cases, the inclusion of different datasets increases the accuracy of prediction. Each dataset must be in a table with rows as samples and features as columns
- Prior knowledge: in order to include prior knowledge in the classification process, the algorithm can include a database that relate the input features in one or several datasets. In this example, we use a pathway database that relates the genes in gene expresion and CNV datasets.
- Output datasets: This algorithm is tunned to 2 class clasification, but is based in SimpleMKL and could be easily modified for regression.

parameters:
- max_unique_kernels: control the max number of features that are included in each kernel.
- nbiter: number of CV iterations
- C :SVM parameter
- kernelt: kernel types that can use the MKL. All the kernel types available are defined in svmkernel.m example: {'gaussian' 'poly'} #gaussian and polinomial kernels
- kerneloptionvect: options of each kernel type. Example = {[0.5 1 2] [1 2]} # gaussian kernels with parameters 0.5, 1 and 2 and polinomial kernel with parameters 1 and 2.
- variablevec: How list of kernels is created. FSMKL options are "ranked" and "rankedoptim". Whole description of all options are in "CreateKernelListWithVariable.m"

note:
For each added datasets, some parameters need to be configured. Follow the examples of PathwaySelectionFSMKL_main.m in order to include new datasets.
The EXP and CNV dataset has prior knowledge (represented by pathways), but clinical data has not. That means that there are one "virtual dataset" for each pathway in CVN and EXP, but just one for each clinical type of data.

output:
The output of the algorithm consist in:
trainR: Train error
testR: Test error
pps: Number of correct matches
ppn: Number of false mathches
auroc: Test AUROC
cm:	test confusion matrix

testRProb: Test error with probability measure
aurocProb: Auroc error with probability measure
cmProb: test confusion matrix with probability measure

testR_2: Test error with outlier removal
auroc_2: Test AUROC with outlier removal
cm_2:	Test confusion matrix with outlier removal

testRProb_2: Test error with outlier removal and probability measure
aurocProb_2: Test AUROC with outlier removal and probability measure
cmProb_2:  Tes confusion matrix with outlier removal and probability measure

Two Pathway files: One with all the samples and other without ourliers, with the selected pathways.

IMPORTANT CONSIDERATIONS:
The size of the final kernel depends of the number of samples and the number of kernels. The number
of kernels dependes of the number of datasets, number of kernel function and number of kernel function options.
N of kernels = 
if rankedoptim: number of dataset * threshold parameter * number of functions * number of options 
if ranked: number of dataset * (threshold parameter-1)! * number of functions * number of options

The size of this kernel can be larger than the memory available in the system becouse the virtualization of the kernel in disk space. 
The number of options should be selected carefully depending of the memory available in the system in order to avoid out of memory errors.


The following filenames are original from SimpleMKL:

- devectorize.c
- devectorize_single.c
- vectorize.c
- vectorize_single.c
- devectorize.dll
- devectorize_single.dll
- vectorize.dll
- vectorize_single.dll
- gradsvmclass.m
- mklkernel.m
- mklsvm.m
- mklsvmupdate.m
- monqp.m
- normalizemeanstd.m
- sumKbeta.m
- svmclass.m
- svmkernel.m
- devectorize.mexa64
- devectorize.mexglx
- devectorize_single.mexa64
- devectorize_single.mexglx
- vectorize.mexa64
- vectorizew.mexglx
- vectorize_single.mexa64
- vectorize_single.mexglx

The following filenames are FS-MKL original or modified from SimpleMKL.
- build_efficientKEFF.m: based in original build_efficientK.m
- CreateKernelListWithVariable.m: based in original
- GetKernelinformationPathway2.m: Original
- getProb.m: Original
- mklkernelEff.m: based in mklkernel.m
- PathwaySelectionFSMKL_main.m
- prepareMultipleDataSet.m
- procesProb.m
- samplingEqual.m
- SaveInformationPathway.m
- UnitTraceNormalization.m: based in original
- GenerateData.m
- pathways.mat
- randomData.mat

The code listed above is distributed under GPL license

The following filenames are original from 
 Cardillo G. (2008) ROC curve: compute a Receiver Operating Characteristics curve.
 http://www.mathworks.com/matlabcentral/fileexchange/19950
 - roc3.m

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 

