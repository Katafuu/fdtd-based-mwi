function available = fbtsParallelAvailable()
%fbtsParallelAvailable A license entitlement alone does not install the toolbox.
available = exist('parpool','file') == 2 && ...
    license('test','Distrib_Computing_Toolbox');
end
