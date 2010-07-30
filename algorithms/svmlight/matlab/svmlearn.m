function model = svmlearn(x,y,parm_string)
% SVMLEARN use SVM-lite to train a Support Vector Machine (SVM)
%   SVMLEARN takes a set of observations (x) as a matrix, and a 
%   a set of labels (y), and a string corresponding to the 
%   command-line arguments for SVM-lite package.
%
%   SVM-lite is compiled as a DLL and is executed through the
%   MEX interface.  There is a penalty incurred through conversion
%   of the vectors from MATLAB's column-major format for arrays to
%   the C standard row-major format.  In practice, this penalty is
%   inconsequential.  The overall performance is nearly identical to
%   that of SVM-Lite running as a stand-alone program.
    clear fun mexsvmclassify;
    clear fun mexsvmlearn;

    try
        model = mexsvmlearn(x,y,parm_string);
    catch
        fprintf(1,'**************************\n');
        lasterror
        fprintf(1,'**************************\n');
        parm_string
        fprintf(1,'**************************\n');
        rethrow(lasterror);
    end
    
    clear fun mexsvmlearn;
    