function Create_Print_Dataset(N,out_csv,seed)

if nargin<1, N=2000; end
if nargin<2, out_csv="jobs_manifest.csv"; end
if nargin<3, seed=123; end
rng(seed);

%% ----------------------- CONFIG ----------------------------------------
materials1 = ["PLA","Tough PLA","ABS","PETG"];
materials2 = ["None","PLA","Tough PLA","ABS","PETG","PVA Natural","Breakaway"];

profiles = ["Balanced","Visual","Draft","Engineering"];
profile_layers = struct();
profile_layers.Balanced    = [0.06 0.10 0.15 0.20 0.30];
profile_layers.Visual      = [0.06 0.10 0.15 0.20];
profile_layers.Draft       = [0.20 0.30];
profile_layers.Engineering = [0.10 0.15 0.20];

infill_patterns   = ["Triangle","Zig zag","gyroid","cubic","grids"];
adhesion_types    = ["none","skirt","brim","raft"];
support_structures= ["tree","normal"];
support_patterns  = ["zigzag","lines","grid","cross","gyroid"];

nozzle_size   = 0.4;
default_speed = 35;

optTemp = containers.Map(["PLA","Tough PLA","ABS","PETG","PVA Natural","Breakaway"],...
                         [200   210         240   240    205          215]);
optBed  = containers.Map(["PLA","Tough PLA","ABS","PETG","PVA Natural","Breakaway"],...
                         [60     60          100   80     45           60 ]);
optFan  = containers.Map(["PLA","Tough PLA","ABS","PETG","PVA Natural","Breakaway"],...
                         [70     60          0     30     60           30 ]);
patternBonus = containers.Map(infill_patterns,[+0.04 +0.02 -0.03 0 -0.01]);

%% ----------------------- TABLE INIT -------------------------------------
T = table();
T.job_id = (1:N).';

T.material_1 = randsample(materials1,N,true).';
T.profile    = randsample(profiles,N,true).';
T.layer_height_mm = zeros(N,1);
for i=1:N
    opts = profile_layers.(char(T.profile(i)));
    T.layer_height_mm(i) = opts(randi(numel(opts)));
end

%% ----------------------- GEOMETRY --------------------------------------
LEN_rng=[20,330]; DEP_rng=[20,240]; HGT_rng=[10,300];
geom_types = ["box","cantilever","U_bracket","spurs","channels"];
geom_overhang_base = containers.Map(geom_types,[0.05 0.45 0.30 0.40 0.25]);

T.model_length_mm=zeros(N,1);
T.model_depth_mm=zeros(N,1);
T.model_height_mm=zeros(N,1);
T.surface_area_mm2=zeros(N,1);
T.volume_mm3=zeros(N,1);
T.compactness=zeros(N,1);
T.aspect_ratio=zeros(N,1);
T.overhang_ratio=zeros(N,1);
T.overhang_area_mm2=zeros(N,1);
T.geometry_type = strings(N,1);
[T.orientation_ex_deg,T.orientation_ey_deg,T.orientation_ez_deg] = deal(zeros(N,1));
T.support_needed = false(N,1);
T.support_structure = strings(N,1);
T.support_pattern   = strings(N,1);
T.support_density_pct = zeros(N,1);

for i=1:N
    L = LEN_rng(1)+rand*(LEN_rng(2)-LEN_rng(1));
    D = DEP_rng(1)+rand*(DEP_rng(2)-DEP_rng(1));
    H = HGT_rng(1)+rand*(HGT_rng(2)-HGT_rng(1));
    T.model_length_mm(i)=L; T.model_depth_mm(i)=D; T.model_height_mm(i)=H;

    SA = 2*(L*D + L*H + D*H);
    V  = L*D*H;
    T.surface_area_mm2(i) = SA;
    T.volume_mm3(i)       = V;
    T.compactness(i)      = V/(SA^1.5);

    dims = sort([L D H]);
    T.aspect_ratio(i) = dims(3)/dims(1);

    gtype = randsample(geom_types,1,true,[.25 .25 .2 .15 .15]);
    T.geometry_type(i) = gtype;

    baseOH = geom_overhang_base(gtype);
    ex = randsample(0:15:165,1); ey = randsample(0:15:165,1); ez = randsample(0:30:330,1);
    T.orientation_ex_deg(i)=ex; T.orientation_ey_deg(i)=ey; T.orientation_ez_deg(i)=ez;

    faceAlignBonus = faceAlignReduction(ex,ey,ez,L,D,H);
    overhang = min(1,max(0, baseOH + 0.12*randn - faceAlignBonus));
    overhang = overhang*(1+0.08*(T.aspect_ratio(i)-1));
    T.overhang_ratio(i)   = overhang;
    T.overhang_area_mm2(i)= overhang*SA;

    T.support_needed(i) = overhang > 0.22;
    if T.support_needed(i)
        if any(gtype==["spurs","channels","cantilever"])
            T.support_structure(i) = "tree";
        else
            T.support_structure(i) = randsample(support_structures,1,true);
        end
        T.support_pattern(i)     = randsample(support_patterns,1,true);
        T.support_density_pct(i) = randi([15 40]);
    else
        T.support_structure(i) = "none";
        T.support_pattern(i)   = "none";
        T.support_density_pct(i)=0;
    end
end

%% ----------------------- MATERIAL 2 (smart) -----------------------------
T.material_2 = strings(N,1);
for i=1:N
    if T.support_needed(i) && rand < 0.7
        % 70% of supported prints use dedicated support material
        if rand < 0.5
            T.material_2(i) = "PVA Natural";
        else
            T.material_2(i) = "Breakaway";
        end
    elseif rand < 0.4
        % 40% of other prints are dual-material (decorative / hybrid)
        T.material_2(i) = randsample(materials2(2:5),1,true); % PLA / Tough PLA / ABS / PETG
    else
        T.material_2(i) = "None";
    end
end

%% ----------------------- MATERIAL COMPATIBILITY -------------------------
compat = string(repmat("Good",N,1));
for i=1:N
    m1=T.material_1(i); m2=T.material_2(i);
    if m2=="None"
        compat(i)="Single Material";
    elseif m1=="ABS" && m2=="PVA Natural"
        compat(i)="Incompatible";
    elseif (m1=="PLA" && (m2=="ABS"||m2=="Tough PLA")) || ...
           (m1=="ABS" && (m2=="PLA"||m2=="Tough PLA"))
        compat(i)="Poor";
    elseif (m1=="PETG" && m2=="PVA Natural") || (m1=="PVA Natural" && m2=="PETG")
        compat(i)="Fair";
    elseif m1==m2
        compat(i)="Excellent";
    else
        compat(i)="Good";
    end
end
T.material_compatibility = compat;

%% ----------------------- PROCESS PARAMS ---------------------------------
[T.nozzle1_temp_c,T.nozzle2_temp_c,T.bed_temp_c] = deal(zeros(N,1));
T.fan_speed_pct = zeros(N,1);
for i=1:N
    m1 = T.material_1(i);
    T.nozzle1_temp_c(i) = round(optTemp(m1)+randn*6);

    if T.material_2(i)=="None"
        T.nozzle2_temp_c(i) = 0;
        T.bed_temp_c(i)     = round(optBed(m1)+randn*3);
    else
        m2 = T.material_2(i);
        T.nozzle2_temp_c(i) = round(optTemp(m2)+randn*6);
        T.bed_temp_c(i)     = round(0.7*optBed(m1)+0.3*optBed(m2)+randn*4);
    end

    if m1=="ABS"
        T.fan_speed_pct(i) = 0;
    else
        T.fan_speed_pct(i) = min(100,max(0, round(optFan(m1)+randn*20)));
    end
end
% Special case: PETG+Breakaway avoid too high fan
for i=1:N
    if (T.material_1(i)=="PETG" && T.material_2(i)=="Breakaway") || ...
       (T.material_2(i)=="PETG" && T.material_1(i)=="Breakaway")
        T.fan_speed_pct(i) = min(T.fan_speed_pct(i),60);
    end
end
T.print_speed_mms = min(120,max(5, round(default_speed + randn(N,1)*12)));

%% ----------------------- INFILL / ADHESION ------------------------------
T.infill_pattern            = randsample(infill_patterns,N,true).';
T.infill_density_pct        = randsample(10:10:100,N,true).';
T.wall_thickness_mm         = nozzle_size*randsample(1:6,N,true).';
T.top_bottom_thickness_mm   = nozzle_size*randsample(1:8,N,true).';
T.adhesion_type             = randsample(adhesion_types,N,true,[0.2 0.25 0.35 0.2]).';

%% ----------------------- FILAMENT CONDITION -----------------------------
T.filament_condition_1 = strings(N,1);
T.filament_condition_2 = strings(N,1);
for i=1:N
    m1=T.material_1(i); m2=T.material_2(i);
    if m1=="PVA Natural"
        T.filament_condition_1(i)=randsample(["good","moist"],1,true,[0.7 0.3]);
    else
        T.filament_condition_1(i)=randsample(["good","aged","moist"],1,true,[0.6 0.25 0.15]);
    end
    if m2=="None"
        T.filament_condition_2(i)="none";
    elseif m2=="PVA Natural"
        T.filament_condition_2(i)=randsample(["good","moist"],1,true,[0.7 0.3]);
    else
        T.filament_condition_2(i)=randsample(["good","aged","moist"],1,true,[0.6 0.25 0.15]);
    end
end
T.prime_tower = (T.material_2 ~= "None");

%% ----------------------- PRINT TIME (hours) -----------------------------
T.print_time_hr = zeros(N,1);
for i=1:N
    baseFill = 0.35*(T.infill_density_pct(i)/100) + 0.15;
    suppFill = 0; if T.support_needed(i), suppFill = 0.10*(T.support_density_pct(i)/30); end
    effFill  = min(0.95, baseFill + suppFill);

    layer_h  = T.layer_height_mm(i);
    flow_rate = 0.4*layer_h*max(10,T.print_speed_mms(i));
    volume_to_print = T.volume_mm3(i)*effFill;

    base_time_hr = (volume_to_print / (flow_rate * 250 * 60)); % hours
    complexity   = 1 + 0.2*(T.support_needed(i)) + 0.3*(T.aspect_ratio(i)/10);
    est_time_hr  = base_time_hr * complexity;

    est_time_hr  = max(1, min(est_time_hr, 120));
    T.print_time_hr(i) = round(est_time_hr * (0.9 + 0.2*rand), 2);
end

%% ----------------------- ORIENTATION CLASS ------------------------------
bins   = [0 15 30 45 60 75 90 360];
labels = {'Flat','LowTilt','MediumTilt','HighTilt','Steep','Upright'};
idx = discretize(T.orientation_ez_deg,bins);
idx(isnan(idx)) = numel(labels);
idx(idx<1 | idx>numel(labels)) = numel(labels);

lbl = repmat(string(labels{end}), height(T),1); % default Upright
for j=1:height(T)
    if idx(j)>=1 && idx(j)<=numel(labels)
        lbl(j) = labels{idx(j)};
    end
end
T.orientation_class = categorical(lbl, labels);

%% ----------------------- WEIGHT (even 0–800 g) --------------------------
rhoFun = @(m) materialDensity(m);
T.weight_mat1_g=zeros(N,1);
T.weight_mat2_g=zeros(N,1);
T.total_weight_g=zeros(N,1);

vol_max = max(T.volume_mm3); if vol_max==0, vol_max=1; end

for i=1:N
    dens1 = rhoFun(T.material_1(i));
    dens2 = 0; if T.material_2(i)~="None", dens2 = rhoFun(T.material_2(i)); end

    f_solid = solidFraction(T.infill_density_pct(i), T.wall_thickness_mm(i), ...
              T.top_bottom_thickness_mm(i), ...
              [T.model_length_mm(i) T.model_depth_mm(i) T.model_height_mm(i)], 0.4);
    V_model = T.volume_mm3(i) * f_solid;

    if T.support_needed(i)
        V1 = 0.9*V_model; V2 = 0.1*V_model*(T.support_density_pct(i)/30);
    else
        V1 = V_model; V2 = 0;
    end

    w1 = dens1*V1; w2 = dens2*V2; raw_total = w1 + w2;

    % Stratified uniform targets to fill 0–800 g evenly
    if i <= round(0.33*N)
        target_wt = 50 + 200*rand;        % 50–250 g
    elseif i <= round(0.66*N)
        target_wt = 250 + 250*rand;       % 250–500 g
    else
        target_wt = 500 + 300*rand;       % 500–800 g
    end

    % Mild correlation to physical volume
    scaleFactor = 0.8 + 0.4 * (T.volume_mm3(i)/vol_max);
    target_wt = min(target_wt * scaleFactor, 800);

    % Preserve material ratios exactly
    if raw_total > 0
        r1 = w1 / raw_total; r2 = w2 / raw_total;
    else
        r1 = 1; r2 = 0;
    end
    w1 = target_wt * r1; w2 = target_wt * r2;

    T.weight_mat1_g(i) = w1;
    T.weight_mat2_g(i) = w2;
    T.total_weight_g(i)= w1 + w2;
end

%% ----------------------- STRENGTH-TO-WEIGHT -----------------------------
T.strength_to_weight = zeros(N,1);
for i=1:N
    geom_boost = 1.0 + 0.1*(T.geometry_type(i)=="box") + 0.05*(T.geometry_type(i)=="U_bracket");
    total_w    = max(1e-6, T.total_weight_g(i));
    strength_factor = (0.35*(T.infill_density_pct(i)/100) + ...
                       0.25*(T.wall_thickness_mm(i)/(3*0.4)) + ...
                       0.25*(T.top_bottom_thickness_mm(i)/(4*0.4))) * geom_boost;
    T.strength_to_weight(i) = strength_factor / total_w;
end

%% ----------------------- SUCCESS MODEL ----------------------------------
compat_bonus = containers.Map(["Excellent","Good","Fair","Poor","Incompatible","Single Material"],...
                              [+0.05       +0.02  0.00  -0.10   -0.25          +0.03]);

T.pass_fail = strings(N,1);
T.print_success_prob = zeros(N,1);
for i=1:N
    p = 0.75;

    ctype = string(T.material_compatibility(i));
    if isKey(compat_bonus,ctype), p = p + compat_bonus(ctype); end

    if T.support_density_pct(i)<10, p=p-0.25;
    elseif T.support_density_pct(i)<15, p=p-0.10; end
    if T.support_density_pct(i)>40, p=p-0.05; end

    p = p + patternBonus(T.infill_pattern(i));

    a = string(T.adhesion_type(i));
    if a=="none", p=p-0.12; elseif a=="skirt", p=p-0.08;
    elseif a=="brim", p=p+0.06; elseif a=="raft", p=p+0.04; end

    fc1=string(T.filament_condition_1(i));
    fc2=string(T.filament_condition_2(i));
    if any([fc1 fc2]=="moist")
        if any([string(T.material_1(i)) string(T.material_2(i))]=="PVA Natural")
            p = p - 0.35;
        else
            p = p - 0.12;
        end
    elseif any([fc1 fc2]=="aged")
        p = p - 0.06;
    end

    density = T.infill_density_pct(i);
    vol_cm3 = T.volume_mm3(i)/1000;
    if     vol_cm3>200 && density<=20, p=p-0.18;
    elseif vol_cm3>100 && density<=20, p=p-0.12;
    elseif vol_cm3>100 && density<=30, p=p-0.06; end

    v = T.print_speed_mms(i);
    if v>=70, p=p-0.15; elseif v>=55, p=p-0.08; elseif v<=15, p=p-0.06; end

    if T.wall_thickness_mm(i)<=0.8,       p=p-0.05; end
    if T.top_bottom_thickness_mm(i)<=0.8, p=p-0.05; end

    if T.support_needed(i) && string(T.support_structure(i))=="tree" && ...
       any(string(T.geometry_type(i))==["spurs","channels","cantilever"])
        p = p + 0.05;
    end

    if T.prime_tower(i), p = p + 0.03; end

    p = min(0.98, max(0.02, p));
    T.print_success_prob(i) = p;
    T.pass_fail(i) = string(ifelse(rand<p,"pass","fail"));
end

%% ----------------------- SAVE CSVs --------------------------------------
% Convert categoricals to strings for CSV
for c=1:width(T)
    if iscategorical(T.(c)), T.(c) = string(T.(c)); end
end

writetable(T, out_csv);
fprintf('Wrote %s with %d rows\n', out_csv, height(T));

% ML subset
features = ["material_1","material_2","material_compatibility","profile","layer_height_mm",...
    "nozzle1_temp_c","nozzle2_temp_c","bed_temp_c","fan_speed_pct","print_speed_mms",...
    "infill_pattern","infill_density_pct","wall_thickness_mm","top_bottom_thickness_mm","adhesion_type",...
    "filament_condition_1","filament_condition_2","prime_tower","geometry_type","model_length_mm","model_depth_mm","model_height_mm",...
    "surface_area_mm2","volume_mm3","compactness","aspect_ratio","overhang_ratio","overhang_area_mm2","support_needed","support_structure",...
    "support_pattern","support_density_pct","orientation_class","print_time_hr","weight_mat1_g","weight_mat2_g","total_weight_g","strength_to_weight",...
    "pass_fail","print_success_prob"];
ML_T = T(:,features);
writetable(ML_T, "print_trainset.csv");
fprintf("Also wrote print_trainset.csv (ML subset)\n");
fprintf("Overall pass rate: %.2f%%\n", 100*mean(ML_T.pass_fail=="pass"));
end

%% ----------------------- HELPERS ---------------------------------------
function bonus = faceAlignReduction(ex,ey,ez,L,D,H)
R = rotz(ez)*roty(ey)*rotx(ex);
zaxis = (R*[0;0;1]).';
[~,upIdx] = max(abs(zaxis));
[~,largestFaceIdx] = max([L*D, L*H, D*H]);
bonus = 0.02;
if (largestFaceIdx==1 && upIdx==3) || (largestFaceIdx~=1 && upIdx~=3)
    bonus = 0.10;
end
end

function M=rotx(a)
c=cosd(a); s=sind(a);
M=[1 0 0; 0 c -s; 0 s c];
end
function M=roty(a)
c=cosd(a); s=sind(a);
M=[c 0 s; 0 1 0; -s 0 c];
end
function M=rotz(a)
c=cosd(a); s=sind(a);
M=[c -s 0; s c 0; 0 0 1];
end

function d=materialDensity(mat)
switch string(mat)
    case "PLA",         d=1.24e-3; % g/mm^3
    case "Tough PLA",   d=1.24e-3;
    case "ABS",         d=1.05e-3;
    case "PETG",        d=1.27e-3;
    case "PVA Natural", d=1.19e-3;
    case "Breakaway",   d=1.20e-3;
    otherwise,          d=1.20e-3;
end
end

function fsol = solidFraction(infill_pct,wall_t,topbot_t,dims,nozzle)
L=dims(1); D=dims(2);
perimeter = 2*(L+D);
wall_lines = max(1, round(wall_t/nozzle));
skin_layers= max(1, round(topbot_t/nozzle));
skin_frac  = min(0.3, 0.02*skin_layers);
wall_frac  = min(0.35, 0.002*wall_lines*perimeter/max(L*D,1));
core_frac  = 0.7*(infill_pct/100);
fsol = min(0.95, max(0.08, core_frac + wall_frac + skin_frac));
end

function y=ifelse(cond,a,b)
if cond, y=a; else, y=b; end
end
