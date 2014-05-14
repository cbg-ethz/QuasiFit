function curateSample(file_name, hamming_cutoff)
switch nargin
    case 0
        fprintf('You need to provide at least a QuasiRecomb FASTA file!\n');
        return;
    case 1
        % default of G_1
        % we assume most users want to use the extended precision sampler
        hamming_cutoff = 1;
end

if hamming_cutoff < 1
    fprintf('k for the hamming cutoff has to be at least 1!\n');
    return;
end

% step 1: load data from QuasiRecomb FASTA
Haplotypes = readGenotypes(file_name);
H = hammingDistance(Haplotypes, hamming_cutoff);

DNA_Haplotypes = cell(1, length(Haplotypes));
for i = 1:length(Haplotypes)
    DNA_Haplotypes{i} = Haplotypes(i).dna;
end
[~, idx_sorted] = sort(DNA_Haplotypes);
Haplotypes = Haplotypes(idx_sorted);

G = sparse(H);
[S, C] = graphconncomp(G);
fprintf('Number of haplotypes in %s: %i\n', file_name, length(Haplotypes));
fprintf('Number of disconnected haplotype clusters; %i\n', S);

if (S == 1)
    fprintf('The haplotype graph G_%i only has one connected component,\njust converting to QuasiFit input file.\n', hamming_cutoff);
end

S_t = S;
i = 1;
while (S_t > 1)
    i = i + 1;
    H_t = hammingDistance(Haplotypes, i);
    G_t = sparse(H_t);
    [S_t, ~] = graphconncomp(G_t);
end
fprintf('k for which G_k has one connected component: %i\n', i);

% step 2: pull connected components apart
CLASSES = {};
for i = 1:S
    CLASSES{i} = Haplotypes(find(C == i));
end

% step 3: span up graph of disconnected components
inserts = [];
testHaplotypes = Haplotypes;

connectedComponentsGraph = zeros(S, S);
connectedComponentsPillars = cell(S, S);

fprintf('Connecting all unconnected components to each other\n');
for i = 1:(S-1)
    for j = (i+1):S
        [DNA_1, DNA_2, MIN] = findBridgePillars(CLASSES{i}, CLASSES{j});
        connectedComponentsGraph(j, i) = MIN;
        connectedComponentsPillars{j, i} = [DNA_1, DNA_2];
    end
end

% step 4: determine minimum spanning tree of graph of disconnected components
[MST, ~] = graphminspantree(sparse(connectedComponentsGraph));
MST = full(MST);

% step 5: insert all unobserved making up the edges of the MST
unobservedHaplotypes = [];
for i = 1:(S-1)
    for j = (i+1):S
        if (MST(j, i) ~= 0)
            tempPillars = connectedComponentsPillars{j, i};
            [constructs] = buildBridge(tempPillars(1), tempPillars(2), hamming_cutoff);
            unobservedHaplotypes = [unobservedHaplotypes constructs];
        end
    end
end

% step 6: sort fullHaplotypes by DNA sequence
fullHaplotypes = [Haplotypes unobservedHaplotypes];
DNA_fullHaplotypes = cell(1, length(fullHaplotypes));

for i = 1:length(fullHaplotypes)
    DNA_fullHaplotypes{i} = fullHaplotypes(i).dna;
end
[~, idx_sorted] = sort(DNA_fullHaplotypes);
fullHaplotypes = fullHaplotypes(idx_sorted);

% step 7: print diagnostic output
H_full = hammingDistance(fullHaplotypes, hamming_cutoff);
[S, C] = graphconncomp(sparse(H_full));

fprintf('=======================================================\n');
fprintf('Number of final haplotypes (including unobserved): %i\n', length(fullHaplotypes));
fprintf('Number of disconnected haplotype clusters: %i\n', S);
fprintf('=======================================================\n');

if (length(fullHaplotypes) > 300)
    fprintf('The final number of haplotypes is more than 300!\nConvergence in the MCMC will be excruciatingly slow,\nkeep this in mind if you proceed with QuasiFit.\n');
end
if (length(unobservedHaplotypes) > 15)
    fprintf('You have more than 15 unobserved haplotypes!\nThis will slow down convergence and decrease\noverall MCMC efficiency.\n');
end

% step 8: write curated output
fprintf('=======================================================\n');
fprintf('Writing output to file\n');
dot_pos = strfind(file_name, '.');
file_stem = file_name(1:(dot_pos-1));

fid=fopen([file_stem, '.in'], 'w');
for i = 1:length(fullHaplotypes)
    fprintf(fid, [fullHaplotypes(i).dna, ',', num2str(fullHaplotypes(i).count), '\n']);
end
fclose=(fid);
end


function [constructs] = buildBridge(DNA_1, DNA_2, hamming_cutoff)
diff_ind = find( DNA_1.dna ~= DNA_2.dna );
constructs = [];
num = length(diff_ind)-1;

for i = hamming_cutoff:hamming_cutoff:num
    temp = DNA_1.dna;
    temp(diff_ind(1:i)) = DNA_2.dna(diff_ind(1:i));
    
    genotype.dna = temp;
    genotype.count = 0;
    constructs = [constructs genotype];
end
end


function [DNA_1, DNA_2, MIN] = findBridgePillars(class_1, class_2)
MIN = 1E7;

for i = 1:length(class_1)
    for j = 1:length(class_2)
        if (sum(class_1(i).dna ~= class_2(j).dna) < MIN)
            DNA_1 = class_1(i);
            DNA_2 = class_2(j);
            MIN = sum(class_1(i).dna ~= class_2(j).dna);
        end
    end
end
end


function [Haplotypes] = readGenotypes(input_file)
f_input = fopen(input_file);
tid = fgetl(f_input);

Haplotypes = [];
while ischar(tid)
    tline = fgetl(f_input);
    freq_pos = strfind(tid, '_');
    
    freq = str2double(tid(freq_pos+1:end));
    
    genotype.dna = tline;
    genotype.count = int32(freq*1E4);
    
    Haplotypes = [Haplotypes genotype];
    
    tid = fgetl(f_input);
end

fclose(f_input);
end


function H = hammingDistance(Genotypes, cutoff)
DIM = length(Genotypes);
H = zeros(DIM, DIM);

for i = 1:DIM
    for j = (i+1):DIM
        dist = sum(Genotypes(i).dna ~= Genotypes(j).dna);
        if (dist > cutoff)
            dist = 0;
        end
        
        H(i,j) = dist;
    end
end

H = H + H';
end