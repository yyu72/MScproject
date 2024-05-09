% Define the data path and file types
datahome = '/home/yi/Documents/MScproject/data';
pathhome = '/home/yi/Documents/MScproject/';

addpath('/home/yi/Documents/MScproject/code');

% Add necessary paths for MATLAB to access additional scripts or functions
addpath(genpath([pathhome, 'Matlab']));  % Adds the Matlab directory and all its subdirectories

% Change the current directory to the data directory
cd(datahome);

% Constants
L = 29903;
file_ext = 'tsv';
filecov_ext = 'bedtools.cov';
prefix = [];
metadataFile = 'metadata.mat';
load(metadataFile);
metadata = T;
clear T;

Patient = unique(metadata.patient_id);
Npatients = numel(Patient);

% filename = cell(1,Npatients);

daymatchstring = 'DAY';
daymatchincrement = 3;

for j = 1:Npatients
    % Loop through each patient by their index from 1 to the total number of patients
    PatientID=Patient{j}
    if ~isempty(prefix)
        % If a filename prefix is specified, execute the following
        filename = [prefix, num2str(Patient{j})];
        % Construct the filename by concatenating the prefix with the patient's ID
        files = dir([filename, '*.', file_ext]);
        % Search for all files matching the constructed filename with the specified file extension (tsv)
        filescov = dir([filename, '*.', filecov_ext]);
        % Similarly, search for all coverage files with the specified coverage file extension (.bedtools.cov)

        for k = 1:numel(files)
            % Iterate over each found file
            indt1 = regexp(files(k).name, daymatchstring) + daymatchincrement;
            % Use a regular expression to find the date in the filename, adjust by daymatchincrement
            indt2 = regexp(files(k).name, file_ext) - 2;
            % Find the position of the file extension and adjust to get the end of the date
            tfile(k) = str2num(files(k).name(indt1:indt2));
            % Extract the date from the filename and convert it to a number, storing it in tfile
        end
    else
        % If no filename prefix is specified, perform the following
        ind = find(strcmp(Patient{j}, metadata.patient_id))
        % Find all entries in the metadata that match the current patient's ID
        fileprefix = metadata.sample_barcode(ind)
        % Retrieve all sample IDs corresponding to the patient

        fs = cell2struct(fileprefix, 'name', 2)
        % Convert the list of sample IDs to a structure array with a field named 'name'
        files = fs;
        % Assign the structure array to files
        filescov = fs;
        % Also assign the structure array to filescov

        skip = false;
        % Initialize a skip flag as false, to determine whether to skip processing for the current patient

        for k = 1:numel(ind)
            % Iterate over all sample IDs matched to the current patient
            files(k).name = [fs(k).name, '.', file_ext];
            % Construct the complete filename by appending the file extension (.tsv)
            filescov(k).name = [fs(k).name, '.', filecov_ext];
            % Construct the coverage file name by appending the coverage file extension (.bedtools.cov)

            if ~isfile(fullfile(datahome, files(k).name))
                % Check if the constructed filename exists in the filesystem
                indname = find(strcmp(fileprefix(k), golbarcode.sample_barcode));
                % If the file does not exist, look up the sample ID in the golbarcode table
                if ~isempty(indname)
                    % If an index was found, perform the following
                    golprefix = golbarcode.gol_number(indname);
                    % Get the corresponding gol number
                    files(k).name = [golprefix, '.', file_ext];
                    % Reconstruct the filename using the gol number
                    filescov(k).name = [golprefix, '.', filecov_ext];
                    % Reconstruct the coverage filename using the gol number
                else
                    % If no valid gol number was found, set the skip flag to true
                    skip = true;
                end
            end
        end

        if ~skip
            % If the skip flag is not set, proceed with the following
            tfile = metadata.date_since_first_pos(ind);
            % Get all dates corresponding to the current patient
            [t, indt] = sort(tfile);
            % Sort the dates, indt is the index after sorting
            files = files(indt);
            % Reorder the files array using the sorted index
            filescov = filescov(indt);
            % Reorder the filescov array using the sorted index
            % Call the function to process the data
            f = NucleotideFrequencies(L, t, files, filescov, PatientID);
        end
    end
end

