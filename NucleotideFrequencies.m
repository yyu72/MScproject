function f=NucleotideFrequencies(L, t, files, filescov, PatientID)

% Sequences contains a cell array of N sequences, where the 1st is the
% reference sequence, to which all other sequences in Sequences have been
% aligned

%This function calculates counts of each type of nucleotide at each
%position in the genome relative to the reference sequences, by reading in
%SNP data from a tsv file from ivar's output 

%This version does not work on vcf files.


tic

%By default the number of nucleotides is presumed to be 6 (4 + ins and deletions) 
nnuc=6;

pvaluethreshold=0.05;

% files = dir([filename,'*.',file_ext])
T = numel(files);%it is presumed that each file corresponds to a different time point
S = cell(1,T);

f = zeros(L,nnuc,T); % frequency array
psequencing = zeros(L,nnuc,T); % pvalues for sequencing error (prob of null model that reads are due to sequencing error)
pstrandbias = ones(L,nnuc,T); % pvalues for strand bias (prob of null model that there is *no* strand bias) 
% — is 0 the correct default?? — have changed to 1 and later REF will be set to NaN

% Qalt = zeros(L,nnuc,T); %Quality score array
% tfile = zeros(1,T);
% daymatchstring = 'DAY'
% daymatchincrement = 3;
% daymatchstring = 'S'
% daymatchincrement = 1;

% filescov = dir([filename,'*.',filescov_ext])
% 
% if numel(filescov)~=T
%     disp('Error number of tsv/vcf files (time-points) files should be the same as number of read-depth/coverage files')
%     return
% end
    
RD = zeros(L,T);


posall = [];
refnucs = cell(1,L);

a=1e5;
% hopen = fopen([filename,'_allelefrequencies.mat']);
% 
% if hopen==-1
    
for k=1:T
    


    if ~strmatch('nofile',files(k).name)

        S{k} = readtable(files(k).name,'FileType','text');
    %         save([filename,'_',num2str(k),'_original_snp_table.mat','Sorig']);
    %         S{k} = Sorig;
    %         [nseg,~] = size(S{k}); %number of segregating position/sites at this timepoint
    
        [nsnps,~] = size(S{k});
        %Add columns for pvalues for sequencing and strandbias
        S{k}.pval_seq = zeros(nsnps,1);
        S{k}.pval_refseq = zeros(nsnps,1);
        S{k}.pval_SB = zeros(nsnps,1);
    
    %     %dir does not list in numerical order so create
    %     indt1=regexp(files(k).name,daymatchstring)+daymatchincrement
    %     indt2=regexp(files(k).name,file_ext)-2
    %     files(k).name(indt1:indt2)
    %     tfile(k) = str2num(files(k).name(indt1:indt2));
    
    %Read-in covergage files with read-depth info
    %     rd = readmatrix(filescov(k).name,'FileType','text');
        rd = readtable(filescov(k).name,'FileType','text');
    %     RD(:,k) = rd(:,3);
        RD(:,k) = rd.Var3;
    
    
        posunique = unique(S{k}.POS); %List of unique sites
        posall = [posall;posunique];
        nseg = numel(posunique)
    
        for n=1:nseg
    %         n
    %         tic
            disp([num2str(n),'/',num2str(nseg)])
            pos = posunique(n); %Find nth unique position
            ind = find(S{k}.POS==pos); %find all entries with this pos 
            % — i.e. at a given position there may be multiple variants and multiple entries
    
    
            Sref = S{k}(ind(1),:); %All n segregating sites at this position should have the same reference allele!
    
            %Calculate totalreads ignoring the total_dp in tsv files, since it
            %doesn't account for ins/dels 
    
            totalreads = Sref.REF_DP;
            for j=1:numel(ind)
                Skj = S{k}(ind(j),:);
    
                totalreads = totalreads + Skj.ALT_DP;
    
            end
    %         totalreads
            %Update the read-depth — N.B. at sites that don't have any SNPs,
            %then the read-depth will be as inputted from cov file—this doesn't
            %always agree with the reads from the tsv files from ivar!
            RD(pos,k) = totalreads;
    
    
    %         count = zeros(1,nnuc);
    %             del=0;
    %             delblock = 0;
    
            %Calculate reference psequencing pvalue
    
            
            refnuc = cell2mat(Sref.REF); %Identity of reference allele
            refnucs{pos} = refnuc;
            refreads = Sref.REF_DP; %Number of reads of the reference allele
            refQ = Sref.REF_QUAL;
            refErr = 10.^(-refQ/10);
    %         totalreads = Sref.TOTAL_DP; %Total number of reads %This is unreliable as total number of reads for ALT +REF doesn't always add up to total_dp!
    
            pseqref = binomialtest(refErr,refreads,totalreads);
    
            if pseqref<-a*eps | pseqref >1+a*eps
                disp('Error REF sequencing p-value not physical')
                save([PatientID,'Error.mat']); 
                return
            end
    
    %         mean_RefSeqErrReads = totalreads*refErr;
    %         
    %         if floor(2*mean_RefSeqErrReads-refreads)>=0
    %             kref1 = 0:floor(2*mean_RefSeqErrReads-refreads):
    %             kref2 = refreads:totalreads;
    %             kref = [kref1,kref2];
    %             pseqref = sum(binopdf(kref,totalreads,refErr));
    %         else
    %             kref = refreads:totalreads;
    %             pseqref = sum(binopdf(kref,totalreads,refErr));
    %         end
    
            switch refnuc
                case {'A','a'}
                    psequencing(pos,1,k) = pseqref;
    %                 psequencing(pos,1,k) = NaN;
                case {'C','c'}
                    psequencing(pos,2,k) =  pseqref;
    %                 psequencing(pos,2,k) = NaN;
                case {'G','g'}
                    psequencing(pos,3,k) =  pseqref;
    %                 psequencing(pos,3,k) = NaN;
                case {'T','t'}
                    psequencing(pos,4,k) =  pseqref;
    %                 psequencing(pos,4,k) = NaN;
                case {'-'}
                    psequencing(pos,5,k) =  pseqref;
    %                 psequencing(pos,5,k) = NaN;
                case {'+'}
                    psequencing(pos,6,k) =  pseqref;
    %                 psequencing(pos,6,k) = NaN;
            end
    
    
            for j=1:numel(ind) %loop over multiple variants whose index is ind
                Skj = S{k}(ind(j),:);
    
                
                altnuc = cell2mat(Skj.ALT); %Identity of alternate allele 
    %                 refnuc = cell2mat(Skj.REF);
                
                altreads = Skj.ALT_DP; %Number of reads of the alternate allele
                
    %                 altfreq = Skj.ALT_FREQ;
    
                
                altQ = Skj.ALT_QUAL;
    %                 refaa = cell2mat(Skj.REF_AA);
    %                 altaa = cell2mat(Skj.ALT_AA);
                refrevreads = Skj.REF_RV;
                altrevreads = Skj.ALT_RV;
    
                
                altErr = 10.^(-altQ/10);
    
                %Filter based on quality score and strand read bias —  
                %Note using a contingency table and Fisher's exact test is 
                %incorrect and too conservative—it is a binomial sampling problem
    
                pseq = binomialtest(altErr,altreads,totalreads);
    
                if pseq<-a*eps | pseq >1+a*eps
                    disp('Error ALT sequencing p-value not physical')
                    save([PatientID,'Error.mat']); 
                    return
                end
    
    %             kalt = altreads:totalreads; %Create vector of reads greater than or equal to observed number altreads
    %             pseq = sum(binopdf(kalt,totalreads,altErr));
                
                %Update table with sequencing pvalue
                S{k}.pval_seq(ind(j))=pseq;
                S{k}.pval_refseq(ind(j))=pseqref;
    
    
                %Use Fisher's exact test for forward vs reverse reads
                %Use a 2-sided test since we are asking symmetrically whether
                %for and reverse read proportions agree, whether one is smaller
                %or greater than the other
                X = [altreads,refreads;altrevreads,refrevreads];
                [~,prevbias,~] = fishertest(X);
    
                if prevbias<-a*eps | prevbias >1+a*eps
                    disp('Error strand-bias p-value not physical')
                    save([PatientID,'Error.mat']); 
                    return
                end
    
                %Update table with strand bias pvalue
                S{k}.pval_SB(ind(j))=prevbias;
    
                switch altnuc(1)
                        case {'A','a'}
                            f(pos,1,k) = altreads/totalreads;
    %                         count(1) = altreads;
                            psequencing(pos,1,k) = pseq;nano
                            pstrandbias(pos,1,k) = prevbias;
    
                        case {'C','c'}
                            f(pos,2,k) = altreads/totalreads;
    %                         count(2) = altreads; 
                            psequencing(pos,2,k) = pseq;
                            pstrandbias(pos,2,k) = prevbias;
                        case {'G','g'}
                            f(pos,3,k) = altreads/totalreads;
    %                         count(3) = altreads;  
                            psequencing(pos,3,k) = pseq;
                            pstrandbias(pos,3,k) = prevbias;
                        case {'T','t'}
                            f(pos,4,k) = altreads/totalreads;
    %                         count(4) = altreads;
                            psequencing(pos,4,k) = pseq;
                            pstrandbias(pos,4,k) = prevbias;
                        case {'-'}
                            f(pos,5,k) = altreads/totalreads;
    %                         count(5) = altreads;
                            psequencing(pos,5,k) = pseq;
                            pstrandbias(pos,5,k) = prevbias;
    %                             del=1;
                        case {'+'}
                            f(pos,6,k) = altreads/totalreads;
                            psequencing(pos,6,k) = pseq;
                            pstrandbias(pos,6,k) = prevbias;
    %                         count(6) = altreads;
                end
    
    
    
            end
    
    
            sumf = sum(f(pos,:,k));
            
            switch refnuc
                case {'A','a'}
                    f(pos,1,k) = 1-sumf;
                case {'C','c'}
                    f(pos,2,k) =  1-sumf;
                case {'G','g'}
                    f(pos,3,k) =  1-sumf;
                case {'T','t'}
                    f(pos,4,k) =  1-sumf;
                case {'-'}
                    f(pos,5,k) =  1-sumf;
                case {'+'}
                    f(pos,6,k) =  1-sumf;
            end
    
    
        end

    else
        t(k) = NaN;
        f(:,:,k) = NaN;
        RD(:,k) = NaN;
        psequencing(:,:,k) = NaN;
        pstrandbias(:,:,k) = NaN;
        S{k} = NaN;


    end

end

%Now at those sites that have segregating variation, but not for all time
%points need to have the correct frequency of the wild type. 
%Loop over all sites and time points again to correct this:

posall =unique(posall); %Unique SNPs across all time points

npos = numel(posall);

for n=1:numel(posall)
    disp([num2str(n),'/',num2str(npos)])

    pos = posall(n);

    for k=1:T

        if ~isnan(t(k))
            x = squeeze(f(pos,:,k));
    
            if x==0
                        
                switch refnucs{pos}
                    case {'A','a'}
                        f(pos,1,k) =  1;
                    case {'C','c'}
                        f(pos,2,k) =  1;
                    case {'G','g'}
                        f(pos,3,k) =  1;
                    case {'T','t'}
                        f(pos,4,k) =  1;
                    case {'-'}
                        f(pos,5,k) =  1;
                    case {'+'}
                        f(pos,6,k) =  1;
                end
    
            end

        end

    end

end

% 
% 
% %Sort the timepoints in ascending order
% [t, indt] = sort(tfile);
% 
% f = f(:,:,indt);
% psequencing = psequencing(:,:,indt);
% pstrandbias = pstrandbias(:,:,indt);
% 
% snptable = S(indt);
% RD = RD(:,indt);

%Remove NaN entries

indt = ~isnan(t);
t= t(indt);
f = f(:,:,indt);
RD = RD(:,indt);
psequencing = psequencing(:,:,indt);
pstrandbias = pstrandbias(:,:,indt);


snptable = S(indt);

mkdir(PatientID)
cd(PatientID)

save([PatientID,'_AlleleFrequencies+snp_table.mat'],'t','f','psequencing','pstrandbias','snptable','RD','refnucs');

cd ..
% end
    
   toc             




end



function p = binomialtest(r,m,ns)

        mean_reads = ns*r;
%         mstar = floor(2*mean_reads-m);

        if mean_reads<m
            mstar = floor(2*mean_reads-m);
            
            if mstar>=0
                k1 = 0:mstar;
                k2 = m:ns;
                k = [k1,k2];
                p = sum(binopdf(k,ns,r));
            elseif mstar==m
                p=1;
            else
                k = m:ns;
                p = sum(binopdf(k,ns,r));
            end

        elseif mean_reads == m
            p=1;
        elseif mean_reads> m
            mstar = ceil(2*mean_reads-m);
            
            if mstar<=ns
                k1 = 0:m;
                k2 = mstar:ns;
                k = [k1,k2];
                p = sum(binopdf(k,ns,r));
            elseif mstar==m
                p=1;
            else
                k = 0:m;
                p = sum(binopdf(k,ns,r));
            end
        end
end










                