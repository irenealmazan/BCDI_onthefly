function [idscns,numscns] = mkidspec(name)
% [idscns,numscns] = mkidspec(name)
% used by readspec to make index for spec file and write to disk
% name = spec data file name (string)
% 13-Aug-96 Brian Stephenson and Carol Thompson ANL
% stephenson@anl.gov 
% elaborated from loadspec.m by Sean Brennan & Anneli Munkholm SSRL

namendx = [name '.ndx']; 
disp(['Creating index file ' namendx]);

fid=fopen(name,'r');
chdat=fread(fid,inf,'uchar');
fclose(fid);

chdat = setstr(chdat)';

% idscns is vector of character indexes to every scan in file

idscns = findstr(chdat,'#S ');

% numscns is vector of scan numbers associated with each scan in file

numscns = zeros(size(idscns));
for ii = 1:length(idscns)
  numscns(ii) = sscanf(chdat(2+idscns(ii):6+idscns(ii)),'%d',1);
end

% Write index to disk

fidi = fopen(namendx,'w');
if fidi < 0
   display('Failed to write index file; no write permission?');
else
   fprintf(fidi,'%i %i\n',[idscns; numscns]);
   fclose(fidi);
end

