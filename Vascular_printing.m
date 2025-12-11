clear all 
% load hollow_lumen_structure.mat or load partitioned_lumen_structure
gfilename = 'Vascular_partitioned_lumen.txt';
fid = fopen(gfilename,'w');
z_height = 0.4 : 0.4 : 0.4 * 21;
for c = 1 : 21
    fprintf(fid,'Z %.3f C %.3f\n',z_height(c),z_height(c));
    shell1 = g_y{c,1};
    shell1 = left_partitioned_lumen_structure{c,1}; % or load left_hollow_lumen_structure
    shell2 = right_partitioned_lumen_structure{c,1}; % or load right_hollow_lumen_structure
    uu1 = length(shell1);
     label = 0;
     long_point(1:1000,1) = 0;
     long_point(1:1000,2) = 0;
uu2 = length(shell2);
if uu1 < uu2
    s = uu1;
    l = uu2;
else
    s = uu2;
    l = uu1;
end
for i = 1 : s
    prez(i,1) = shell1(i,1);
    prez(i,2) = shell1(i,2);
    prey(i,1) = shell2(i,1);
    prey(i,2) = shell2(i,2);
end
dnan1 = isnan(prez);
dnan2 = isnan(prey);
jnan1 = find(dnan1(:,1)>0);
jnan2 = find(dnan2(:,1)>0);
postz = prez;
posty = prey;
juu1 = 0;
juu2 = 0;
for i = 1 : length(jnan1)
    prez(jnan1(i)-juu1,:) = [];
    juu1 = juu1+1;
end
for i = 1 : length(jnan2)
    prey(jnan2(i)-juu2,:) = [];
    juu2 = juu2+1;
end
z_new = prez;
y_new = prey;
z = length(prez);
y = length(prey);
if z > y
    s = y;
    l = z;
else
    s = z;
    l = y;
end
for i = 1 : s
    listw(i,1) = z_new(i,1);
    listw(i,2) = z_new(i,2);
    listw(i,3) = y_new(i,1);
    listw(i,4) = y_new(i,2);
end
if isempty(jnan1)
    jnan1(1) = 0;
end
if isempty(jnan2)
    jnan2(1) = 0;
end
if jnan1(1) > jnan2(1)
    s = jnan2(1);
    l = jnan1(1);
else
    s = jnan1(1);
    l = jnan2(1);
end
for i  = 1 : s
    listy(i,1) = shell1(i,1);
    listy(i,2) = shell1(i,2);
    listy(i,3) = shell2(i,1);
    listy(i,4) = shell2(i,2);
end
if length(z_new) > length(y_new)
    s = y_new;
    l = z_new;
else
    s = z_new;
    l = y_new;
end
tep = zeros(length(s),1);
list_tep = zeros(length(s),4);
time_tep = 0;
for i = 1 : length(s)-1
    distancez = ((z_new(i+1,1)-z_new(i,1))^2+(z_new(i+1,2)-z_new(i,2))^2)^0.5;
    distancey = ((y_new(i+1,1)-y_new(i,1))^2+(y_new(i+1,2)-y_new(i,2))^2)^0.5;
    ddlist(i,1) = distancez;
    ddlist(i,2) = distancey;
    if abs(((distancey-distancez) / 1)) > 2 
        time_tep = time_tep + 1;
        tep(time_tep,1) = i;
        list_tep(time_tep,1) = z_new(i,1);
        list_tep(time_tep,2) = z_new(i,2);
        list_tep(time_tep,3) = y_new(i,1);
        list_tep(time_tep,4) = y_new(i,2);
    end
end
a = 0;
for i = 1 : size(list_tep,1) * size(list_tep,2)
    if list_tep(i) ~= 0
        a = a + 1;
    end
end
if a ~= 0
for i = 1 : time_tep
    jl(i,1) = list_tep(i,1);
    jl(i,2) = list_tep(i,2);
    jl(i,3) = list_tep(i,3); 
    jl(i,4) = list_tep(i,4);
end
for i =  1 : length(shell1)
    for j = 1 : size(jl,1)
        if jl(j,1) == shell1(i,1) && jl(j,2) == shell1(i,2)
            long_point(j,1) = i;
        end
    end
end
for i = 1 : length(shell2)
    for j = 1 : size(jl,1)
        if jl(j,3) == shell2(i,1) && jl(j,4) == shell2(i,2)
            long_point(j,2) = i;
        end
    end
end
for i = 1 : time_tep
    jl(i,1) = list_tep(i,1); 
    jl(i,2) = list_tep(i,2);
    jl(i,3) = list_tep(i,3); 
    jl(i,4) = list_tep(i,4);
end
for i =  1 : length(shell1)
    for j = 1 : size(jl,1)
        if jl(j,1) == shell1(i,1) && jl(j,2) == shell1(i,2)
            long_point(j,1) = i;
        end
    end
end
for i = 1 : length(shell2)
    for j = 1 : size(jl,1)
        if jl(j,3) == shell2(i,1) && jl(j,4) == shell2(i,2)
            long_point(j,2) = i;
        end
    end
end
for i =1 : size(long_point,1)-1 
    if long_point(i+1,1) < long_point(i,1)
        for i1 = 1 : length(jnan1) 
            juli1(i,i1) = long_point(i,1)-jnan1(i1,1);
            juli1_abs(i,i1) = abs(long_point(i,1)-jnan1(i1,1));
            min_juli1 = min(juli1_abs); 
            [row,col] = find(juli1_abs==min_juli1);
            if long_point(i,1) <= max(jnan1) 
                if col == 1
                    long_point(i,1) = 1;
                end
            elseif long_point(i,1) > max(jnan1) 
                long_point(i,1) = max(jnan1)+1;
            end
        end
    elseif long_point(i+1,2) < long_point(i,2)
        for i2 = 1 : length(jnan2)
            juli2(i,i2) = long_point(i,2)-jnan2(i2,1);
            juli2_abs(i,i2) = abs(long_point(i,2)-jnan2(i2,1));
            min_juli2 = min(juli2_abs);
            [row,col] = find(juli2_abs==min_juli2);
            if long_point(i,2) <= max(jnan2)
                if col == 1
                    long_point(i,2) = 1;
                else
                    long_point(i,2) = jnan2(col,1)-jnan2(col-1,1)+1;
                end
            elseif long_point(i,1) > max(jnan2)
                long_point(i,2) = max(jnan2) + 1;
            end
        end
    end
end
bianchang = size(long_point,1);
for i = 1 : bianchang
    lp1(i,1) = long_point(i,1);
    lp2(i,1) = long_point(i,2);
end
lp = [lp1;lp2];
lp_tep = unique(lp);
long_point = lp_tep;
if size(long_point,1) == 1
    label = 1;
end
long_point(find(long_point==0)) = [];
unique(long_point);
if isempty(long_point)
    long_point(1,1) = 0;
end
shelly = shell1;
shellz = shell2;
fprintf(fid,'X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',shell1(1,1)-shell2(1,1),shell1(1,2),z_height(c),-shell2(1,1),shell2(1,2),z_height(c));
u1 = length(shelly); 
u2 = length(shellz);
du1 = 0;
du2 = 0;
uu = 0;
uu_bc = 0;
overlap = zeros(50);
if u1 > u2
    s = shellz;
    t = shelly;
else
    s = shelly;
    t = shellz;
end
flag = 1;
AAA = size(s,1);
long_point(find(long_point==AAA)) = [];
unique(long_point);
for i_point = 1 + uu_bc : size(s,1)
    zsx = shellz(i_point+du1,1);
    zsy = shellz(i_point+du1,2);
    ysx = shelly(i_point+du2,1);
    ysy = shelly(i_point+du2,2);
    du1 = 0;
    du2 = 0;
    if ~ismember(i_point,long_point(:,1)) 
        if ~isnan(zsx) && ~isnan(ysx)
            if flag
              if i_point ~= 2
                fprintf(fid,'X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                flag = 0;
              else
                 fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                 flag = 0;
              end
            end
        if flag == 0
            if i_point == 2
                fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                X1 = ysx-zsx;
                X2 = ysy;
                A1 = -zsx;
                A2 = zsy;
            else
                fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                X1 = ysx-zsx;
                X2 = ysy;
                A1 = -zsx;
                A2 = zsy;
                NEVER = 1;
            end
        else
            continue
        end
    elseif ~isnan(zsx) && isnan(ysx)
        fprintf(fid,'G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',shelly(i_point+1,1)-zsx,shelly(i_point+1,2),z_height(c),-zsx,zsy,z_height(c));
        X1 = shelly(i_point+1,1)-zsx;
        X2 = shelly(i_point+1,2);
        A1 = -zsx;
        A2 = zsy;
    elseif isnan(zsx) && ~isnan(ysx) 
        fprintf(fid,'G1 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-shellz(i_point+1,1),ysy,z_height(c),-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
        X1 = ysx-shellz(i_point+1,1);
        X2 = ysy;
        A1 = -shellz(i_point+1,1);
        A2 = shellz(i_point+1,2);
        elseif isnan(zsx) && isnan(ysx)
            fprintf(fid,'X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',shelly(i_point+1,1)-shellz(i_point+1,1),shelly(i_point+1,2),z_height(c),-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
            X1 = shelly(i_point+1,1)-shellz(i_point+1,1);
            X2 = shelly(i_point+1,2);
            A1 = -shellz(i_point+1,1);
            A2 = shellz(i_point+1,2);
        end
    elseif ismember(i_point,long_point(:,1)) 
        if i_point == long_point(1,1) && label == 0 
            if ~ismember(i_point+1,long_point(:,1)) 
                if ~isnan(shelly(i_point+1,1)) && ~isnan(shellz(i_point+1,1))
                    fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                    fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
                    fprintf(fid,'X %.3f Y %.3f Z %.3f\n',(ysx-zsx)+(-shellz(i_point+1,1)+zsx),ysy,z_height(c));
                    fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',(ysx-zsx)+(-shellz(i_point+1,1)+zsx),ysy,z_height(c));                    
                    fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',shelly(i_point+1,1)-shellz(i_point+1,1),shelly(i_point+1,2),z_height(c));
                    X1 = shelly(i_point+1,1)-shellz(i_point+1,1);
                    X2 = shelly(i_point+1,2);
                    A1 = -shellz(i_point+1,1);
                    A2 = shellz(i_point+1,2);
                elseif isnan(shelly(i_point+1,1)) && ~isnan(shellz(i_point+1,1)) 
                    fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                    fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
                    fprintf(fid,'X %.3f Y %.3f Z %.3f\n',((ysx-zsx)+(shelly(i_point+2,1)-shellz(i_point+1,1)))*0.5,(ysy+shelly(i_point+2,2))*0.5,z_height(c));
                    fprintf(fid,'X %.3f Y %.3f Z %.3f\n',shelly(i_point+2,1)-shellz(i_point+1,1),shelly(i_point+2,2),z_height(c));
                    X1 = shelly(i_point+2,1)-shellz(i_point+1,1);
                    X2 = shelly(i_point+2,2);
                    A1 = -shellz(i_point+1,1);
                    A2 = shellz(i_point+1,2);
                elseif ~isnan(shelly(i_point+1,1)) && isnan(shellz(i_point+1,1)) 
                    fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                    fprintf(fid,'A %.3f B %.3f C %.3f\n',-shellz(i_point+2,1),shellz(i_point+2,2),z_height(c));
                    fprintf(fid,'X %.3f Y %.3f Z %3.f\n',(ysx-zsx)+(-shellz(i_point+2,1)+zsx),ysy,z_height(c));
                    fprintf(fid,'G1 X %.3f Y %.3f Z %3.f\n',(ysx-zsx)+(-shellz(i_point+2,1)+zsx),ysy,z_height(c));
                    fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',shelly(i_point+1,1)-shellz(i_point+2,1),shelly(i_point+1,2),z_height(c));
                    X1 = shelly(i_point+1,1)-shellz(i_point+2,1);
                    X2 = shelly(i_point+1,2);
                    A1 = -shellz(i_point+2,1);
                    A2 = shellz(i_point+2,2);
                end
            elseif ismember(i_point+1,long_point(:,1)) 
                for cycle = 1 : 1000
                    if ismember(i_point+cycle,long_point(:,1))
                        uu = uu + 1; 
                        uu_bc = uu_bc + 1;
                    else
                        break
                    end
                end
                if ~isnan(ysx) && ~isnan(zsx)
                fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                X1 = ysx-zsx;
                X2 = ysy;
                A1 = -zsx;
                A2 = zsy;
                elseif ~isnan(ysx) && isnan(zsx)
                    fprintf(fid,'G1 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-shellz(i_point+1,1),ysy,z_height(c),-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
                    X1 = ysx-shellz(i_point+1,1);
                    X2 = ysy;
                    A1 = -shellz(i_point+1,1);
                    A2 = shellz(i_point+1,2);
                elseif isnan(ysx) && ~isnan(zsx)
                    fprintf(fid,'G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',shelly(i_point+1,1)-zsx,shelly(i_point+1,2),z_height(c),-zsx,zsy,z_height(c));
                    X1 = shelly(i_point+1,1)-zsx;
                    X2 = shelly(i_point+1,2);
                    A1 = -zsx;
                    A2 = zsy;
                end
                for i = 1 : uu
                    if ~isnan(shelly(i_point+i,1)) && ~isnan(shellz(i_point+i,1))
                        if ~isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i,1)+shellz(i_point+i-1,1)),shelly(i_point+i-1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i,1)+shellz(i_point+i-1,1)),shelly(i_point+i-1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',shelly(i_point+i,1)-shellz(i_point+i,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        elseif isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',((shelly(i_point+i-2,1)-shellz(i_point+i-1,1))+(shelly(i_point+i,1)-shellz(i_point+i,1)))*0.5,(shelly(i_point+i-2,2)+shelly(i_point+i,2))*0.5,z_height(c));
                            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',shelly(i_point+i,1)-shellz(i_point+i,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        elseif ~isnan(shelly(i_point+i-1,1)) && isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f B %.3f C %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-2,1))+(-shellz(i_point+i,1)+shellz(i_point+i-2,1)),shelly(i_point+i-1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-2,1))+(-shellz(i_point+i,1)+shellz(i_point+i-2,1)),shelly(i_point+i-1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',shelly(i_point+i,1)-shellz(i_point+i,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        elseif isnan(shelly(i_point+i-1,1)) && isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f B %.3f C %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',shelly(i_point+i,1)-shellz(i_point+i,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        end
                    elseif ~isnan(shelly(i_point+i,1)) && isnan(shellz(i_point+i,1)) 
                        if ~isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f B %.3f C %.3f\n',-shellz(i_point+i+1,1),shellz(i_point+i+1,2),z_height(c));
                            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i+1,1)+shellz(i_point+i-1,1)),shelly(i_point+i-1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i+1,1)+shellz(i_point+i-1,1)),shelly(i_point+i-1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',shelly(i_point+i,1)-shellz(i_point+i+1,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i+1,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i+1,1);
                            A2 = shellz(i_point+i+1,2);
                            overlap(i,1) = i_point + i;
                        elseif isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f B %.3f C %.3f\n',-shellz(i_point+i+1,1),shellz(i_point+i+1,2),z_height(c));
                            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',(shelly(i_point+i-2,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i+1,1)+shellz(i_point+i-1,1)),shelly(i_point+i-2,2),z_height(c));
                            fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',(shelly(i_point+i-2,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i+1,1)+shellz(i_point+i-1,1)),shelly(i_point+i-2,2),z_height(c));
                            fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',shelly(i_point+i,1)-shellz(i_point+i+1,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1) - shellz(i_point+i+1,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i+1,1);
                            A2 = shellz(i_point+i+1,2);
                            overlap(i,1) = i_point + i;
                        end
                    elseif isnan(shelly(i_point+i,1)) && ~isnan(shellz(i_point+i,1))
                        if ~isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',shelly(i_point+i+1,1)-shellz(i_point+i,1),shelly(i_point+i+1,2),z_height(c));
                            X1 = shelly(i_point+i+1,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i+1,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        elseif ~isnan(shelly(i_point+i-1,1)) && isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',shelly(i_point+i+1,1)-shellz(i_point+i,1),shelly(i_point+i+1,2),z_height(c));
                            X1 = shelly(i_point+i+1,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i+1,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        end
                    elseif isnan(shelly(i_point+i,1)) && isnan(shellz(i_point+i,1))
                        if ~isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f B %.3f C %.3f\n',-shellz(i_point+i+1,1),shellz(i_point+i+1,2),z_height(c));
                            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',shelly(i_point+i+1,1)-shellz(i_point+i+1,1),shelly(i_point+i+1,2),z_height(c));
                            X1 = shelly(i_point+i+1,1)-shellz(i_point+i+1,1);
                            X2 = shelly(i_point+i+1,2);
                            A1 = -shellz(i_point+i+1,1);
                            A2 = shellz(i_point+i+1,2);
                            overlap(i,1) = i_point + i;
                        end
                    end
                end
                uu = 0;
            end
        end 
if i_point ~= long_point(end,1) && i_point ~= long_point(1,1) 
        if ~ismember(i_point+1,long_point(:,1)) && ~ismember(i_point,overlap(:,1))
            if ~isnan(ysx) && ~isnan(zsx)
                H1 = ysx-zsx;
                H11 = ysy;
                H2 = -zsx;
                H22 = zsy;
            elseif isnan(ysx) && ~isnan(zsx)
                H1 = shelly(i_point+1,1)-zsx;
                H11 = shelly(i_point+1,2);
                H2 = -zsx;
                H22 = zsy;
            elseif ~isnan(ysx) && isnan(zsx)
                H1 = ysx-shellz(i_point+1,1);
                H11 = ysy;
                H2 = -shellz(i_point+1,1);
                H22 = shellz(i_point+1,2);
            end
                if ~isnan(shelly(i_point+1,1)) && ~isnan(shellz(i_point+1,1))
                    fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',H1,H11,z_height(c),H2,H22,z_height(c));
                    fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
                    fprintf(fid,'X %.3f Y %.3f Z %.3f\n',H1+(-shellz(i_point+1,1)-H2),H11,z_height(c));
                    fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',H1+(-shellz(i_point+1,1)-H2),H11,z_height(c));
                    fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',shelly(i_point+1,1)-shellz(i_point+1,1),shelly(i_point+1,2),z_height(c));
                    X1 = shelly(i_point+1,1)-shellz(i_point+1,1);
                    X2 = shelly(i_point+1,2);
                    A1 = -shellz(i_point+1,1);
                    A2 = shellz(i_point+1,2);
                elseif isnan(shelly(i_point+1,1)) && ~isnan(shellz(i_point+1,1)) 
                    fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',H1,H11,z_height(c),H2,H22,z_height(c));
                    fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
                    fprintf(fid,'X %.3f Y %.3f Z %.3f\n',((H1)+(shelly(i_point+2,1)-shellz(i_point+1,1)))*0.5,(H11+shelly(i_point+2,2))*0.5,z_height(c));
                    fprintf(fid,'X %.3f Y %.3f Z %.3f\n',shelly(i_point+2,1)-shellz(i_point+1,1),shelly(i_point+2,2),z_height(c));
                    X1 = shelly(i_point+2,1)-shellz(i_point+1,1);
                    X2 = shelly(i_point+2,2);
                    A1 = -shellz(i_point+1,1);
                    A2 = shellz(i_point+1,2);
                elseif ~isnan(shelly(i_point+1,1)) && isnan(shellz(i_point+1,1)) 
                    fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',H1,H11,z_height(c),H2,H22,z_height(c));
                    fprintf(fid,'A %.3f B %.3f C %.3f\n',-shellz(i_point+2,1),shellz(i_point+2,2),z_height(c));
                    fprintf(fid,'X %.3f Y %.3f Z %.3f\n',H1+(-shellz(i_point+2,1)-H2),H11,z_height(c));
                    fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',H1+(-shellz(i_point+2,1)-H2),H11,z_height(c));
                    fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',shelly(i_point+1,1)-shellz(i_point+2,1),shelly(i_point+1,2),z_height(c));
                    X1 = shelly(i_point+1,1)-shellz(i_point+2,1);
                    X2 = shelly(i_point+1,2);
                    A1 = -shellz(i_point+2,1);
                    A2 = shellz(i_point+2,2);
                end
        end 
        if ismember(i_point+1,long_point(:,1)) 
            if ismember(i_point,overlap(:,1)) 
                never = 1;
            elseif ~ismember(i_point,overlap(:,1))
                for cycle = 1 : 1000
                    if ismember(i_point+cycle,long_point(:,1))
                        uu = uu + 1;
                        uu_bc = uu_bc + 1;
                    else
                        break
                    end
                end 
               if ~isnan(ysx) && ~isnan(zsx)
                fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                X1 = ysx-zsx;
                X2 = ysy;
                A1 = -zsx;
                A2 = zsy;
                elseif ~isnan(ysx) && isnan(zsx)
                    fprintf(fid,'G1 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-shellz(i_point+1,1),ysy,z_height(c),-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
                    X1 = ysx-shellz(i_point+1,1);
                    X2 = ysy;
                    A1 = -shellz(i_point+1,1);
                    A2 = shellz(i_point+1,2);
                elseif isnan(ysx) && ~isnan(zsx)
                    fprintf(fid,'G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',shelly(i_point+1,1)-zsx,shelly(i_point+1,2),z_height(c),-zsx,zsy,z_height(c));
                    X1 = shelly(i_point+1,1)-zsx;
                    X2 = shelly(i_point+1,2);
                    A1 = -zsx;
                    A2 = zsy;
                elseif isnan(ysx) && isnan(zsx)
                    fprintf(fid,'X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',shelly(i_point+1,1)-shellz(i_point+1,1),shelly(i_point+1,2),z_height(c),-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
                    X1 = shelly(i_point+1,1)-shellz(i_point+1,1);
                    X2 = shelly(i_point+1,2);
                    A1 = -shellz(i_point+1,1);
                    A2 = shellz(i_point+1,2);
                end
                for i = 1 : uu
                    if ~isnan(shelly(i_point+i,1)) && ~isnan(shellz(i_point+i,1)) 
                        if ~isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                        fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                        fprintf(fid,'X %.3f Y %.3f Z %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i,1)+shellz(i_point+i-1,1)),shelly(i_point+i-1,2),z_height(c));
                        fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i,1)+shellz(i_point+i-1,1)),shelly(i_point+i-1,2),z_height(c));
                        fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',shelly(i_point+i,1)-shellz(i_point+i,1),shelly(i_point+i,2),z_height(c));
                        X1 = shelly(i_point+i,1)-shellz(i_point+i,1);
                        X2 = shelly(i_point+i,2);
                        A1 = -shellz(i_point+i,1);
                        A2 = shellz(i_point+i,2);
                        overlap(i,1) = i_point + i;
                        elseif isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',((shelly(i_point+i-2,1)-shellz(i_point+i-1,1))+(shelly(i_point+i,1)-shellz(i_point+i,1)))*0.5,(shelly(i_point+i-2,2)+shelly(i_point+i,2))*0.5,z_height(c));
                            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',shelly(i_point+i,1)-shellz(i_point+i,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        elseif ~isnan(shelly(i_point+i-1,1)) && isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f B %.3f C %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-2,1))+(-shellz(i_point+i,1)+shellz(i_point+i-2,1)),shelly(i_point+i-1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-2,1))+(-shellz(i_point+i,1)+shellz(i_point+i-2,1)),shelly(i_point+i-1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',shelly(i_point+i,1)-shellz(i_point+i,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        elseif isnan(shelly(i_point+i-1,1)) && isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f B %.3f C %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',shelly(i_point+i,1)-shellz(i_point+i,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        end
                    elseif ~isnan(shelly(i_point+i,1)) && isnan(shellz(i_point+i,1)) 
                        if ~isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                        fprintf(fid,'A %.3f B %.3f C %.3f\n',-shellz(i_point+i+1,1),shellz(i_point+i+1,2),z_height(c));
                        fprintf(fid,'X %.3f Y %.3f Z %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i+1,1)+shellz(i_point+i-1,1)),shelly(i_point+i-1,2),z_height(c));
                        fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i+1,1)+shellz(i_point+i-1,1)),shelly(i_point+i-1,2),z_height(c));
                        fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',shelly(i_point+i,1)-shellz(i_point+i+1,1),shelly(i_point+i,2),z_height(c));
                        X1 = shelly(i_point+i,1)-shellz(i_point+i+1,1);
                        X2 = shelly(i_point+i,2);
                        A1 = -shellz(i_point+i+1,1);
                        A2 = shellz(i_point+i+1,2);
                        overlap(i,1) = i_point + i;
                        elseif isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f B %.3f C %.3f\n',-shellz(i_point+i+1,1),shellz(i_point+i+1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',shelly(i_point+i,1)-shellz(i_point+i+1,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i+1,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i+1,1);
                            A2 = shellz(i_point+i+1,2);
                            overlap(i,1) = i_point + i;
                        end
                    elseif isnan(shelly(i_point+i,1)) && ~isnan(shellz(i_point+i,1))
                         if ~isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',shelly(i_point+i+1,1)-shellz(i_point+i,1),shelly(i_point+i+1,2),z_height(c));
                            X1 = shelly(i_point+i+1,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i+1,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        elseif ~isnan(shelly(i_point+i-1,1)) && isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',shelly(i_point+i+1,1)-shellz(i_point+i,1),shelly(i_point+i+1,2),z_height(c));
                            X1 = shelly(i_point+i+1,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i+1,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        end
                    elseif isnan(shelly(i_point+i,1)) && isnan(shellz(i_point+i,1))
                            if ~isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f B %.3f C %.3f\n',-shellz(i_point+i+1,1),shellz(i_point+i+1,2),z_height(c));
                            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',shelly(i_point+i+1,1)-shellz(i_point+i+1,1),shelly(i_point+i+1,2),z_height(c));
                            X1 = shelly(i_point+i+1,1)-shellz(i_point+i+1,1);
                            X2 = shelly(i_point+i+1,2);
                            A1 = -shellz(i_point+i+1,1);
                            A2 = shellz(i_point+i+1,2);
                            overlap(i,1) = i_point + i;
                            end
                    end
                end
                uu = 0;
            end
        end 
end
if i_point == long_point(end,1) && (label == 1 || label == 0) 
     if ~isnan(shelly(i_point,1)) && ~isnan(shellz(i_point,1))
        fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
        H1 = ysx-zsx;
        H11 = ysy;
        H2 = -zsx;
        H22 = zsy;
    elseif isnan(shelly(i_point,1)) && ~isnan(shellz(i_point,1))
        fprintf(fid,'G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',shelly(i_point+1,1)-zsx,shelly(i_point+1,2),z_height(c),-zsx,zsy,z_height(c));
        H1 = shelly(i_point+1,1)-zsx;
        H11 = shelly(i_point+1,2);
        H2 = -zsx;
        H22 = zsy;
    elseif ~isnan(shelly(i_point,1)) && isnan(shellz(i_point,1))
        fprintf(fid,'G1 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-shellz(i_point+1,1),ysy,z_height(c),-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
        H1 = ysx-shellz(i_point+1,1);
        H11 = ysy;
        H2 = -shellz(i_point+1,1);
        H22 = shellz(i_point+1,2);
    end
    if ~ismember(i_point-1,long_point(:,1)) 
        if ~isnan(shelly(i_point+1,1)) && ~isnan(shellz(i_point+1,1)) 
            fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',H1,H11,z_height(c),H2,H22,z_height(c));
            fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',H1+(-shellz(i_point+1,1)-H2),H11,z_height(c));
            fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',H1+(-shellz(i_point+1,1)-H2),H11,z_height(c));
            fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',shelly(i_point+1,1)-shellz(i_point+1,1),shelly(i_point+1,2),z_height(c));
            X1 = shelly(i_point+1,1)-shellz(i_point+1,1);
            X2 = shelly(i_point+1,2);
            A1 = -shellz(i_point+1,1);
            A2 = shellz(i_point+1,2);
        elseif isnan(shelly(i_point+1,1)) && ~isnan(shellz(i_point+1,1))
            fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',H1,H11,z_height(c),H2,H22,z_height(c));
            fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',((H1)+(shelly(i_point+2,1)-shellz(i_point+1,1)))*0.5,(H11+shelly(i_point+2,2))*0.5,z_height(c));
            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',shelly(i_point+2,1)-shellz(i_point+1,1),shelly(i_point+2,2),z_height(c));
            X1 = shelly(i_point+2,1)-shellz(i_point+1,1);
            X2 = shelly(i_point+2,2);
            A1 = -shellz(i_point+1,1);
            A2 = shellz(i_point+1,2);
        elseif ~isnan(shelly(i_point+1,1)) && isnan(shellz(i_point+1,1)) 
            fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',H1,H11,z_height(c),H2,H22,z_height(c));
            fprintf(fid,'A %.3f B %.3f C %.3f\n',-shellz(i_point+2,1),shellz(i_point+2,2),z_height(c));
            fprintf(fid,'X %.3f Y %.3f Z %.3f\n',H1+(-shellz(i_point+2,1)-H2),H11,z_height(c));
            fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',H1+(-shellz(i_point+2,1)-H2),H11,z_height(c));
            fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',shelly(i_point+1,1)-shellz(i_point+2,1),shelly(i_point+1,2),z_height(c));
            X1 = shelly(i_point+1,1)-shellz(i_point+2,1);
            X2 = shelly(i_point+1,2);
            A1 = -shellz(i_point+2,1);
            A2 = shellz(i_point+2,2);
        end
    elseif ismember(i_point-1,long_point(:,1)) 
        never = 1;
    end 
end
    end
end
fprintf(fid,'X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',X1,-20,50,A1,-20,50);
if u1 ~= u2
    flag = 1;
    for i_point = size(s,1) : size(t,1)
        pointx = t(i_point,1);
        pointy = t(i_point,2);
        if u1 < u2 
            if ~isnan(pointx)
                if flag
                    if i_point ~= size(s,1)
                        fprintf(fid,'A %.3f B %.3f C %.3f\n',A1,A2,z_height(c));
                        fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',A1,A2,z_height(c));
                        flag = 0;
                    end
                end
                if flag == 0
                    if i_point == size(s,1)
                        fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-pointx,pointy,z_height(c));
                        A1 = -pointx;
                    else
                        if ~isnan(pointx)
                        fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-pointx,pointy,z_height(c));
                        A1 = -pointx;
                        nnnn = 1;
                        end
                    end
                else
                    continue
                end
            else 
                fprintf(fid,'A %.3f B %.3f C %.3f\n',-t(i_point+1,1),t(i_point+1,2),z_height(c));
                A1 = -t(i_point+1,1);
            end
        else 
            if ~isnan(pointx)
                if flag
                    if i_point ~= size(s,1)
                        fprintf(fid,'X %.3f Y %.3f Z %.3f\n',X1,X2,z_height(c));
                        fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',X1,X2,z_height(c));
                        flag = 0;
                    end
                end
                if flag == 0
                    if i_point == size(s,1)
                        fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',pointx+A1,pointy,z_height(c));
                        X1 = pointx + A1;
                    else
                        fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',pointx+A1,pointy,z_height(c));
                        X1 = pointx + A1;
                        nnnn = 0;
                    end
                else 
                    continue
                end
            else 
                fprintf(fid,'X %.3f Y %.3f Z %.3f\n',t(i_point+1,1)+A1,t(i_point+1,2),z_height(c));
                X1 = t(i_point+1,1)+A1;
            end
        end
    end
    if nnnn == 1
        fprintf(fid,'A %.3f B %.3f C %.3f\n',A1,-20,50);
    elseif nnnn == 0
        fprintf(fid,'X %.3f Y %.3f Z %.3f\n',X1,-20,50);
    end
else
    time = 1;
end
elseif a == 0 
shelly = shell1;
shellz = shell2;
fprintf(fid,'X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',shelly(1,1)-shellz(1,1),shelly(1,2),z_height(c),-shellz(1,1),shellz(1,2),z_height(c));
u1 = length(shelly);
u2 = length(shellz);
if u1 > u2
    s = shellz;
    t = shelly;
else
    s = shelly;
    t = shellz;
end
flag = 1;
for i_point = 2 : size(s,1)
    zsx = shellz(i_point,1);
    zsy = shellz(i_point,2);
    ysx = shelly(i_point,1);
    ysy = shelly(i_point,2);
    if ~isnan(zsx) && ~isnan(ysx)
        if flag
            if i_point ~= 2
                fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                flag = 0;
            end
        end
        if flag == 0
            if i_point == 2
                fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                X1 = ysx-zsx;
                X2 = ysy;
                A1 = -zsx;
                A2 = zsy;
            else
                fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                X1 = ysx-zsx;
                X2 = ysy;
                A1 = -zsx;
                A2 = zsy;
                never = -1;
            end
        else
            continue
        end
    elseif ~isnan(zsx) && isnan(ysx)
        fprintf(fid,'G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',shelly(i_point+1,1)-zsx,shelly(i_point+1,2),z_height(c),-zsx,zsy,z_height(c));
        X1 = shelly(i_point+1,1)-zsx;
        X2 = shelly(i_point+1,2);
        A1 = -zsx;
        A2 = zsy;
    elseif isnan(zsx) && ~isnan(ysx)
        fprintf(fid,'G1 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',ysx-shellz(i_point+1,1),ysy,z_height(c),-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
        X1 = ysx-shellz(i_point+1,1);
        X2 = ysy;
        A1 = -shellz(i_point+1,1);
        A2 = shellz(i_point+1,2);
    else
        fprintf(fid,'X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',shelly(i_point+1,1)-shellz(i_point+1,1),shelly(i_point+1,2),z_height(c),-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
    end
end
fprintf(fid,'X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',X1,-20,50,A1,-20,50);
if u1 ~= u2
    flag = 1;
    for i_point = size(s,1) : size(t,1)
        pointx = t(i_point,1);
        pointy = t(i_point,2);
        if u1 < u2 
            if ~isnan(pointx)
                if flag
                    if i_point ~= size(s,1)
                        fprintf(fid,'A %.3f B %.3f C %.3f\n',A1,A2,z_height(c));
                        fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',A1,A2,z_height(c));
                        flag = 0;
                    end
                end
                if flag == 0
                    if i_point == size(s,1)
                        fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-pointx,pointy,z_height(c));
                        A1 = -pointx;
                    else
                        if ~isnan(pointx)
                        fprintf(fid,'G2 A %.3f B %.3f C %.3f\n',-pointx,pointy,z_height(c));
                        A1 = -pointx;
                        nnnn = 1;
                        end
                    end
                else
                    continue
                end
            else 
                fprintf(fid,'A %.3f B %.3f C %.3f\n',-t(i_point+1,1),t(i_point+1,2),z_height(c));
                A1 = -t(i_point+1,1);
            end
        else 
            if ~isnan(pointx)
                if flag
                    if i_point ~= size(s,1)
                        fprintf(fid,'X %.3f Y %.3f Z %.3f\n',X1,X2,z_height(c));
                        fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',X1,X2,z_height(c));
                        flag = 0;
                    end
                end
                if flag == 0
                    if i_point == size(s,1)
                        fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',pointx+A1,pointy,z_height(c));
                        X1 = pointx + A1;
                    else
                        fprintf(fid,'G1 X %.3f Y %.3f Z %.3f\n',pointx+A1,pointy,z_height(c));
                        X1 = pointx + A1;
                        nnnn = 0;
                    end
                else 
                    continue
                end
            else 
                fprintf(fid,'X %.3f Y %.3f Z %.3f\n',t(i_point+1,1)+A1,t(i_point+1,2),z_height(c));
                X1 = t(i_point+1,1)+A1;
            end
        end
    end
    if nnnn == 1
        fprintf(fid,'A %.3f B %.3f C %.3f\n',A1,-20,50);
    elseif nnnn == 0
        fprintf(fid,'X %.3f Y %.3f Z %.3f\n',X1,-20,50);
    end
else
    time = 1;
end
end
%% 
    fprintf(fid,'T %.3f E %.3f\n',z_height(c),z_height(c));
    shell1 = left_partitioned_lumen_sacrificial_structure{c,1}; % or load left_hollow_lumen_sacrificial_structure
    shell2 = right_partitioned_lumen_sacrificial_structure{c,1}; % or load right_hollow_lumen_sacrificial_structure
    uu1 = length(shell1);
     label = 0;
     long_point(1:1000,1) = 0;
     long_point(1:1000,2) = 0;
uu2 = length(shell2);
if uu1 < uu2
    s = uu1;
    l = uu2;
else
    s = uu2;
    l = uu1;
end
for i = 1 : s
    prez(i,1) = shell1(i,1);
    prez(i,2) = shell1(i,2);
    prey(i,1) = shell2(i,1);
    prey(i,2) = shell2(i,2);
end
dnan1 = isnan(prez);
dnan2 = isnan(prey);
jnan1 = find(dnan1(:,1)>0);
jnan2 = find(dnan2(:,1)>0);
postz = prez;
posty = prey;
juu1 = 0;
juu2 = 0;
for i = 1 : length(jnan1)
    prez(jnan1(i)-juu1,:) = [];
    juu1 = juu1+1;
end
for i = 1 : length(jnan2)
    prey(jnan2(i)-juu2,:) = [];
    juu2 = juu2+1;
end
z_new = prez;
y_new = prey;
z = length(prez);
y = length(prey);
if z > y
    s = y;
    l = z;
else
    s = z;
    l = y;
end
for i = 1 : s
    listw(i,1) = z_new(i,1);
    listw(i,2) = z_new(i,2);
    listw(i,3) = y_new(i,1);
    listw(i,4) = y_new(i,2);
end
if isempty(jnan1)
    jnan1(1) = 0;
end
if isempty(jnan2)
    jnan2(1) = 0;
end
if jnan1(1) > jnan2(1)
    s = jnan2(1);
    l = jnan1(1);
else
    s = jnan1(1);
    l = jnan2(1);
end
for i  = 1 : s
    listy(i,1) = shell1(i,1);
    listy(i,2) = shell1(i,2);
    listy(i,3) = shell2(i,1);
    listy(i,4) = shell2(i,2);
end
if length(z_new) > length(y_new)
    s = y_new;
    l = z_new;
else
    s = z_new;
    l = y_new;
end
tep = zeros(length(s),1);
list_tep = zeros(length(s),4);
time_tep = 0;
for i = 1 : length(s)-1
    distancez = ((z_new(i+1,1)-z_new(i,1))^2+(z_new(i+1,2)-z_new(i,2))^2)^0.5;
    distancey = ((y_new(i+1,1)-y_new(i,1))^2+(y_new(i+1,2)-y_new(i,2))^2)^0.5;
    ddlist(i,1) = distancez;
    ddlist(i,2) = distancey;
    if abs(((distancey-distancez) / 1)) > 2 
        time_tep = time_tep + 1;
        tep(time_tep,1) = i;
        list_tep(time_tep,1) = z_new(i,1);
        list_tep(time_tep,2) = z_new(i,2);
        list_tep(time_tep,3) = y_new(i,1);
        list_tep(time_tep,4) = y_new(i,2);
    end
end
a = 0;
for i = 1 : size(list_tep,1) * size(list_tep,2)
    if list_tep(i) ~= 0
        a = a + 1;
    end
end
if a ~= 0
for i = 1 : time_tep
    jl(i,1) = list_tep(i,1); 
    jl(i,2) = list_tep(i,2);
    jl(i,3) = list_tep(i,3); 
    jl(i,4) = list_tep(i,4);
end
for i =  1 : length(shell1)
    for j = 1 : size(jl,1)
        if jl(j,1) == shell1(i,1) && jl(j,2) == shell1(i,2)
            long_point(j,1) = i;
        end
    end
end
for i = 1 : length(shell2)
    for j = 1 : size(jl,1)
        if jl(j,3) == shell2(i,1) && jl(j,4) == shell2(i,2)
            long_point(j,2) = i;
        end
    end
end
for i = 1 : time_tep
    jl(i,1) = list_tep(i,1); 
    jl(i,2) = list_tep(i,2);
    jl(i,3) = list_tep(i,3); 
    jl(i,4) = list_tep(i,4);
end
for i =  1 : length(shell1)
    for j = 1 : size(jl,1)
        if jl(j,1) == shell1(i,1) && jl(j,2) == shell1(i,2)
            long_point(j,1) = i;
        end
    end
end
for i = 1 : length(shell2)
    for j = 1 : size(jl,1)
        if jl(j,3) == shell2(i,1) && jl(j,4) == shell2(i,2)
            long_point(j,2) = i;
        end
    end
end
for i =1 : size(long_point,1)-1 
    if long_point(i+1,1) < long_point(i,1)
        for i1 = 1 : length(jnan1) 
            juli1(i,i1) = long_point(i,1)-jnan1(i1,1);
            juli1_abs(i,i1) = abs(long_point(i,1)-jnan1(i1,1));
            min_juli1 = min(juli1_abs); 
            [row,col] = find(juli1_abs==min_juli1);
            if long_point(i,1) <= max(jnan1) 
                if col == 1
                    long_point(i,1) = 1;
                end
            elseif long_point(i,1) > max(jnan1)
                long_point(i,1) = max(jnan1)+1;
            end
        end
    elseif long_point(i+1,2) < long_point(i,2)
        for i2 = 1 : length(jnan2)
            juli2(i,i2) = long_point(i,2)-jnan2(i2,1);
            juli2_abs(i,i2) = abs(long_point(i,2)-jnan2(i2,1));
            min_juli2 = min(juli2_abs);
            [row,col] = find(juli2_abs==min_juli2);
            if long_point(i,2) <= max(jnan2)
                if col == 1
                    long_point(i,2) = 1;
                else
                    long_point(i,2) = jnan2(col,1)-jnan2(col-1,1)+1;
                end
            elseif long_point(i,1) > max(jnan2)
                long_point(i,2) = max(jnan2) + 1;
            end
        end
    end
end
bianchang = size(long_point,1);
for i = 1 : bianchang
    lp1(i,1) = long_point(i,1);
    lp2(i,1) = long_point(i,2);
end
lp = [lp1;lp2];
lp_tep = unique(lp);
long_point = lp_tep;
if size(long_point,1) == 1
    label = 1;
end
long_point(find(long_point==0)) = [];
unique(long_point);
if isempty(long_point)
    long_point(1,1) = 0;
end
shelly = shell1;
shellz = shell2;
fprintf(fid,'X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',shell1(1,1)-shell2(1,1),shell1(1,2),z_height(c),-shell2(1,1),shell2(1,2),z_height(c));
u1 = length(shelly); 
u2 = length(shellz);
du1 = 0;
du2 = 0;
uu = 0;
uu_bc = 0;
overlap = zeros(50);
if u1 > u2 
    s = shellz;
    t = shelly;
else
    s = shelly;
    t = shellz;
end
flag = 1;
AAA = size(s,1);
long_point(find(long_point==AAA)) = [];
unique(long_point);
for i_point = 1 + uu_bc : size(s,1)
    zsx = shellz(i_point+du1,1);
    zsy = shellz(i_point+du1,2);
    ysx = shelly(i_point+du2,1);
    ysy = shelly(i_point+du2,2);
    du1 = 0;
    du2 = 0;
    if ~ismember(i_point,long_point(:,1)) 
        if ~isnan(zsx) && ~isnan(ysx)
            if flag
              if i_point ~= 2
                fprintf(fid,'X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                flag = 0;
              else
                 fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                 flag = 0;
              end
            end
        if flag == 0
            if i_point == 2
                fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                X1 = ysx-zsx;
                X2 = ysy;
                A1 = -zsx;
                A2 = zsy;
            else
                fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                X1 = ysx-zsx;
                X2 = ysy;
                A1 = -zsx;
                A2 = zsy;
                NEVER = 1;
            end
        else
            continue
        end
    elseif ~isnan(zsx) && isnan(ysx) 
        fprintf(fid,'G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',shelly(i_point+1,1)-zsx,shelly(i_point+1,2),z_height(c),-zsx,zsy,z_height(c));
        X1 = shelly(i_point+1,1)-zsx;
        X2 = shelly(i_point+1,2);
        A1 = -zsx;
        A2 = zsy;
    elseif isnan(zsx) && ~isnan(ysx) 
        fprintf(fid,'G1 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-shellz(i_point+1,1),ysy,z_height(c),-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
        X1 = ysx-shellz(i_point+1,1);
        X2 = ysy;
        A1 = -shellz(i_point+1,1);
        A2 = shellz(i_point+1,2);
        elseif isnan(zsx) && isnan(ysx)
            fprintf(fid,'X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',shelly(i_point+1,1)-shellz(i_point+1,1),shelly(i_point+1,2),z_height(c),-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
            X1 = shelly(i_point+1,1)-shellz(i_point+1,1);
            X2 = shelly(i_point+1,2);
            A1 = -shellz(i_point+1,1);
            A2 = shellz(i_point+1,2);
        end
    elseif ismember(i_point,long_point(:,1)) 
        if i_point == long_point(1,1) && label == 0 
            if ~ismember(i_point+1,long_point(:,1)) 
                if ~isnan(shelly(i_point+1,1)) && ~isnan(shellz(i_point+1,1))
                    fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                    fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
                    fprintf(fid,'X %.3f S %.3f T %.3f\n',(ysx-zsx)+(-shellz(i_point+1,1)+zsx),ysy,z_height(c));
                    fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',(ysx-zsx)+(-shellz(i_point+1,1)+zsx),ysy,z_height(c));                    
                    fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',shelly(i_point+1,1)-shellz(i_point+1,1),shelly(i_point+1,2),z_height(c));
                    X1 = shelly(i_point+1,1)-shellz(i_point+1,1);
                    X2 = shelly(i_point+1,2);
                    A1 = -shellz(i_point+1,1);
                    A2 = shellz(i_point+1,2);
                elseif isnan(shelly(i_point+1,1)) && ~isnan(shellz(i_point+1,1)) 
                    fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                    fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
                    fprintf(fid,'X %.3f S %.3f T %.3f\n',((ysx-zsx)+(shelly(i_point+2,1)-shellz(i_point+1,1)))*0.5,(ysy+shelly(i_point+2,2))*0.5,z_height(c));
                    fprintf(fid,'X %.3f S %.3f T %.3f\n',shelly(i_point+2,1)-shellz(i_point+1,1),shelly(i_point+2,2),z_height(c));
                    X1 = shelly(i_point+2,1)-shellz(i_point+1,1);
                    X2 = shelly(i_point+2,2);
                    A1 = -shellz(i_point+1,1);
                    A2 = shellz(i_point+1,2);
                elseif ~isnan(shelly(i_point+1,1)) && isnan(shellz(i_point+1,1)) 
                    fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                    fprintf(fid,'A %.3f D %.3f E %.3f\n',-shellz(i_point+2,1),shellz(i_point+2,2),z_height(c));
                    fprintf(fid,'X %.3f S %.3f T %3.f\n',(ysx-zsx)+(-shellz(i_point+2,1)+zsx),ysy,z_height(c));
                    fprintf(fid,'G1 X %.3f S %.3f T %3.f\n',(ysx-zsx)+(-shellz(i_point+2,1)+zsx),ysy,z_height(c));
                    fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',shelly(i_point+1,1)-shellz(i_point+2,1),shelly(i_point+1,2),z_height(c));
                    X1 = shelly(i_point+1,1)-shellz(i_point+2,1);
                    X2 = shelly(i_point+1,2);
                    A1 = -shellz(i_point+2,1);
                    A2 = shellz(i_point+2,2);
                end
            elseif ismember(i_point+1,long_point(:,1)) 
                for cycle = 1 : 1000
                    if ismember(i_point+cycle,long_point(:,1))
                        uu = uu + 1; 
                        uu_bc = uu_bc + 1;
                    else
                        break
                    end
                end
                if ~isnan(ysx) && ~isnan(zsx)
                fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                X1 = ysx-zsx;
                X2 = ysy;
                A1 = -zsx;
                A2 = zsy;
                elseif ~isnan(ysx) && isnan(zsx)
                    fprintf(fid,'G1 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-shellz(i_point+1,1),ysy,z_height(c),-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
                    X1 = ysx-shellz(i_point+1,1);
                    X2 = ysy;
                    A1 = -shellz(i_point+1,1);
                    A2 = shellz(i_point+1,2);
                elseif isnan(ysx) && ~isnan(zsx)
                    fprintf(fid,'G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',shelly(i_point+1,1)-zsx,shelly(i_point+1,2),z_height(c),-zsx,zsy,z_height(c));
                    X1 = shelly(i_point+1,1)-zsx;
                    X2 = shelly(i_point+1,2);
                    A1 = -zsx;
                    A2 = zsy;
                end
                for i = 1 : uu
                    if ~isnan(shelly(i_point+i,1)) && ~isnan(shellz(i_point+i,1))
                        if ~isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f S %.3f T %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i,1)+shellz(i_point+i-1,1)),shelly(i_point+i-1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i,1)+shellz(i_point+i-1,1)),shelly(i_point+i-1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',shelly(i_point+i,1)-shellz(i_point+i,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        elseif isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f S %.3f T %.3f\n',((shelly(i_point+i-2,1)-shellz(i_point+i-1,1))+(shelly(i_point+i,1)-shellz(i_point+i,1)))*0.5,(shelly(i_point+i-2,2)+shelly(i_point+i,2))*0.5,z_height(c));
                            fprintf(fid,'X %.3f S %.3f T %.3f\n',shelly(i_point+i,1)-shellz(i_point+i,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        elseif ~isnan(shelly(i_point+i-1,1)) && isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f D %.3f E %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f S %.3f T %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-2,1))+(-shellz(i_point+i,1)+shellz(i_point+i-2,1)),shelly(i_point+i-1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-2,1))+(-shellz(i_point+i,1)+shellz(i_point+i-2,1)),shelly(i_point+i-1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',shelly(i_point+i,1)-shellz(i_point+i,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        elseif isnan(shelly(i_point+i-1,1)) && isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f D %.3f E %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f S %.3f T %.3f\n',shelly(i_point+i,1)-shellz(i_point+i,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        end
                    elseif ~isnan(shelly(i_point+i,1)) && isnan(shellz(i_point+i,1)) 
                        if ~isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f D %.3f E %.3f\n',-shellz(i_point+i+1,1),shellz(i_point+i+1,2),z_height(c));
                            fprintf(fid,'X %.3f S %.3f T %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i+1,1)+shellz(i_point+i-1,1)),shelly(i_point+i-1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i+1,1)+shellz(i_point+i-1,1)),shelly(i_point+i-1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',shelly(i_point+i,1)-shellz(i_point+i+1,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i+1,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i+1,1);
                            A2 = shellz(i_point+i+1,2);
                            overlap(i,1) = i_point + i;
                        elseif isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f D %.3f E %.3f\n',-shellz(i_point+i+1,1),shellz(i_point+i+1,2),z_height(c));
                            fprintf(fid,'X %.3f S %.3f T %.3f\n',(shelly(i_point+i-2,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i+1,1)+shellz(i_point+i-1,1)),shelly(i_point+i-2,2),z_height(c));
                            fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',(shelly(i_point+i-2,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i+1,1)+shellz(i_point+i-1,1)),shelly(i_point+i-2,2),z_height(c));
                            fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',shelly(i_point+i,1)-shellz(i_point+i+1,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1) - shellz(i_point+i+1,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i+1,1);
                            A2 = shellz(i_point+i+1,2);
                            overlap(i,1) = i_point + i;
                        end
                    elseif isnan(shelly(i_point+i,1)) && ~isnan(shellz(i_point+i,1))
                        if ~isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f S %.3f T %.3f\n',shelly(i_point+i+1,1)-shellz(i_point+i,1),shelly(i_point+i+1,2),z_height(c));
                            X1 = shelly(i_point+i+1,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i+1,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        elseif ~isnan(shelly(i_point+i-1,1)) && isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f S %.3f T %.3f\n',shelly(i_point+i+1,1)-shellz(i_point+i,1),shelly(i_point+i+1,2),z_height(c));
                            X1 = shelly(i_point+i+1,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i+1,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        end
                    elseif isnan(shelly(i_point+i,1)) && isnan(shellz(i_point+i,1))
                        if ~isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f D %.3f E %.3f\n',-shellz(i_point+i+1,1),shellz(i_point+i+1,2),z_height(c));
                            fprintf(fid,'X %.3f S %.3f T %.3f\n',shelly(i_point+i+1,1)-shellz(i_point+i+1,1),shelly(i_point+i+1,2),z_height(c));
                            X1 = shelly(i_point+i+1,1)-shellz(i_point+i+1,1);
                            X2 = shelly(i_point+i+1,2);
                            A1 = -shellz(i_point+i+1,1);
                            A2 = shellz(i_point+i+1,2);
                            overlap(i,1) = i_point + i;
                        end
                    end
                end
                uu = 0;
            end
        end 
if i_point ~= long_point(end,1) && i_point ~= long_point(1,1) 
        if ~ismember(i_point+1,long_point(:,1)) && ~ismember(i_point,overlap(:,1))  
            if ~isnan(ysx) && ~isnan(zsx)
                H1 = ysx-zsx;
                H11 = ysy;
                H2 = -zsx;
                H22 = zsy;
            elseif isnan(ysx) && ~isnan(zsx)
                H1 = shelly(i_point+1,1)-zsx;
                H11 = shelly(i_point+1,2);
                H2 = -zsx;
                H22 = zsy;
            elseif ~isnan(ysx) && isnan(zsx)
                H1 = ysx-shellz(i_point+1,1);
                H11 = ysy;
                H2 = -shellz(i_point+1,1);
                H22 = shellz(i_point+1,2);
            end
                if ~isnan(shelly(i_point+1,1)) && ~isnan(shellz(i_point+1,1))
                    fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',H1,H11,z_height(c),H2,H22,z_height(c));
                    fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
                    fprintf(fid,'X %.3f S %.3f T %.3f\n',H1+(-shellz(i_point+1,1)-H2),H11,z_height(c));
                    fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',H1+(-shellz(i_point+1,1)-H2),H11,z_height(c));
                    fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',shelly(i_point+1,1)-shellz(i_point+1,1),shelly(i_point+1,2),z_height(c));
                    X1 = shelly(i_point+1,1)-shellz(i_point+1,1);
                    X2 = shelly(i_point+1,2);
                    A1 = -shellz(i_point+1,1);
                    A2 = shellz(i_point+1,2);
                elseif isnan(shelly(i_point+1,1)) && ~isnan(shellz(i_point+1,1)) 
                    fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',H1,H11,z_height(c),H2,H22,z_height(c));
                    fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
                    fprintf(fid,'X %.3f S %.3f T %.3f\n',((H1)+(shelly(i_point+2,1)-shellz(i_point+1,1)))*0.5,(H11+shelly(i_point+2,2))*0.5,z_height(c));
                    fprintf(fid,'X %.3f S %.3f T %.3f\n',shelly(i_point+2,1)-shellz(i_point+1,1),shelly(i_point+2,2),z_height(c));
                    X1 = shelly(i_point+2,1)-shellz(i_point+1,1);
                    X2 = shelly(i_point+2,2);
                    A1 = -shellz(i_point+1,1);
                    A2 = shellz(i_point+1,2);
                elseif ~isnan(shelly(i_point+1,1)) && isnan(shellz(i_point+1,1)) 
                    fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',H1,H11,z_height(c),H2,H22,z_height(c));
                    fprintf(fid,'A %.3f D %.3f E %.3f\n',-shellz(i_point+2,1),shellz(i_point+2,2),z_height(c));
                    fprintf(fid,'X %.3f S %.3f T %.3f\n',H1+(-shellz(i_point+2,1)-H2),H11,z_height(c));
                    fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',H1+(-shellz(i_point+2,1)-H2),H11,z_height(c));
                    fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',shelly(i_point+1,1)-shellz(i_point+2,1),shelly(i_point+1,2),z_height(c));
                    X1 = shelly(i_point+1,1)-shellz(i_point+2,1);
                    X2 = shelly(i_point+1,2);
                    A1 = -shellz(i_point+2,1);
                    A2 = shellz(i_point+2,2);
                end
        end 
        if ismember(i_point+1,long_point(:,1)) 
            if ismember(i_point,overlap(:,1)) 
                never = 1;
            elseif ~ismember(i_point,overlap(:,1))
                for cycle = 1 : 1000
                    if ismember(i_point+cycle,long_point(:,1))
                        uu = uu + 1;
                        uu_bc = uu_bc + 1;
                    else
                        break
                    end
                end 
               if ~isnan(ysx) && ~isnan(zsx)
                fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                X1 = ysx-zsx;
                X2 = ysy;
                A1 = -zsx;
                A2 = zsy;
                elseif ~isnan(ysx) && isnan(zsx)
                    fprintf(fid,'G1 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-shellz(i_point+1,1),ysy,z_height(c),-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
                    X1 = ysx-shellz(i_point+1,1);
                    X2 = ysy;
                    A1 = -shellz(i_point+1,1);
                    A2 = shellz(i_point+1,2);
                elseif isnan(ysx) && ~isnan(zsx)
                    fprintf(fid,'G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',shelly(i_point+1,1)-zsx,shelly(i_point+1,2),z_height(c),-zsx,zsy,z_height(c));
                    X1 = shelly(i_point+1,1)-zsx;
                    X2 = shelly(i_point+1,2);
                    A1 = -zsx;
                    A2 = zsy;
                elseif isnan(ysx) && isnan(zsx)
                    fprintf(fid,'X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',shelly(i_point+1,1)-shellz(i_point+1,1),shelly(i_point+1,2),z_height(c),-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
                    X1 = shelly(i_point+1,1)-shellz(i_point+1,1);
                    X2 = shelly(i_point+1,2);
                    A1 = -shellz(i_point+1,1);
                    A2 = shellz(i_point+1,2);
                end
                for i = 1 : uu
                    if ~isnan(shelly(i_point+i,1)) && ~isnan(shellz(i_point+i,1)) 
                        if ~isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                        fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                        fprintf(fid,'X %.3f S %.3f T %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i,1)+shellz(i_point+i-1,1)),shelly(i_point+i-1,2),z_height(c));
                        fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i,1)+shellz(i_point+i-1,1)),shelly(i_point+i-1,2),z_height(c));
                        fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',shelly(i_point+i,1)-shellz(i_point+i,1),shelly(i_point+i,2),z_height(c));
                        X1 = shelly(i_point+i,1)-shellz(i_point+i,1);
                        X2 = shelly(i_point+i,2);
                        A1 = -shellz(i_point+i,1);
                        A2 = shellz(i_point+i,2);
                        overlap(i,1) = i_point + i;
                        elseif isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f S %.3f T %.3f\n',((shelly(i_point+i-2,1)-shellz(i_point+i-1,1))+(shelly(i_point+i,1)-shellz(i_point+i,1)))*0.5,(shelly(i_point+i-2,2)+shelly(i_point+i,2))*0.5,z_height(c));
                            fprintf(fid,'X %.3f S %.3f T %.3f\n',shelly(i_point+i,1)-shellz(i_point+i,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        elseif ~isnan(shelly(i_point+i-1,1)) && isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f D %.3f E %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f S %.3f T %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-2,1))+(-shellz(i_point+i,1)+shellz(i_point+i-2,1)),shelly(i_point+i-1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-2,1))+(-shellz(i_point+i,1)+shellz(i_point+i-2,1)),shelly(i_point+i-1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',shelly(i_point+i,1)-shellz(i_point+i,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        elseif isnan(shelly(i_point+i-1,1)) && isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f D %.3f E %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f S %.3f T %.3f\n',shelly(i_point+i,1)-shellz(i_point+i,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        end
                    elseif ~isnan(shelly(i_point+i,1)) && isnan(shellz(i_point+i,1)) 
                        if ~isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                        fprintf(fid,'A %.3f D %.3f E %.3f\n',-shellz(i_point+i+1,1),shellz(i_point+i+1,2),z_height(c));
                        fprintf(fid,'X %.3f S %.3f T %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i+1,1)+shellz(i_point+i-1,1)),shelly(i_point+i-1,2),z_height(c));
                        fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',(shelly(i_point+i-1,1)-shellz(i_point+i-1,1))+(-shellz(i_point+i+1,1)+shellz(i_point+i-1,1)),shelly(i_point+i-1,2),z_height(c));
                        fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',shelly(i_point+i,1)-shellz(i_point+i+1,1),shelly(i_point+i,2),z_height(c));
                        X1 = shelly(i_point+i,1)-shellz(i_point+i+1,1);
                        X2 = shelly(i_point+i,2);
                        A1 = -shellz(i_point+i+1,1);
                        A2 = shellz(i_point+i+1,2);
                        overlap(i,1) = i_point + i;
                        elseif isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f D %.3f E %.3f\n',-shellz(i_point+i+1,1),shellz(i_point+i+1,2),z_height(c));
                            fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',shelly(i_point+i,1)-shellz(i_point+i+1,1),shelly(i_point+i,2),z_height(c));
                            X1 = shelly(i_point+i,1)-shellz(i_point+i+1,1);
                            X2 = shelly(i_point+i,2);
                            A1 = -shellz(i_point+i+1,1);
                            A2 = shellz(i_point+i+1,2);
                            overlap(i,1) = i_point + i;
                        end
                    elseif isnan(shelly(i_point+i,1)) && ~isnan(shellz(i_point+i,1))
                         if ~isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f S %.3f T %.3f\n',shelly(i_point+i+1,1)-shellz(i_point+i,1),shelly(i_point+i+1,2),z_height(c));
                            X1 = shelly(i_point+i+1,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i+1,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        elseif ~isnan(shelly(i_point+i-1,1)) && isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-shellz(i_point+i,1),shellz(i_point+i,2),z_height(c));
                            fprintf(fid,'X %.3f S %.3f T %.3f\n',shelly(i_point+i+1,1)-shellz(i_point+i,1),shelly(i_point+i+1,2),z_height(c));
                            X1 = shelly(i_point+i+1,1)-shellz(i_point+i,1);
                            X2 = shelly(i_point+i+1,2);
                            A1 = -shellz(i_point+i,1);
                            A2 = shellz(i_point+i,2);
                            overlap(i,1) = i_point + i;
                        end
                    elseif isnan(shelly(i_point+i,1)) && isnan(shellz(i_point+i,1))
                            if ~isnan(shelly(i_point+i-1,1)) && ~isnan(shellz(i_point+i-1,1))
                            fprintf(fid,'A %.3f D %.3f E %.3f\n',-shellz(i_point+i+1,1),shellz(i_point+i+1,2),z_height(c));
                            fprintf(fid,'X %.3f S %.3f T %.3f\n',shelly(i_point+i+1,1)-shellz(i_point+i+1,1),shelly(i_point+i+1,2),z_height(c));
                            X1 = shelly(i_point+i+1,1)-shellz(i_point+i+1,1);
                            X2 = shelly(i_point+i+1,2);
                            A1 = -shellz(i_point+i+1,1);
                            A2 = shellz(i_point+i+1,2);
                            overlap(i,1) = i_point + i;
                            end
                    end
                end
                uu = 0;
            end
        end 
end 
if i_point == long_point(end,1) && (label == 1 || label == 0) 
     if ~isnan(shelly(i_point,1)) && ~isnan(shellz(i_point,1))
        fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
        H1 = ysx-zsx;
        H11 = ysy;
        H2 = -zsx;
        H22 = zsy;
    elseif isnan(shelly(i_point,1)) && ~isnan(shellz(i_point,1))
        fprintf(fid,'G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',shelly(i_point+1,1)-zsx,shelly(i_point+1,2),z_height(c),-zsx,zsy,z_height(c));
        H1 = shelly(i_point+1,1)-zsx;
        H11 = shelly(i_point+1,2);
        H2 = -zsx;
        H22 = zsy;
    elseif ~isnan(shelly(i_point,1)) && isnan(shellz(i_point,1))
        fprintf(fid,'G1 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-shellz(i_point+1,1),ysy,z_height(c),-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
        H1 = ysx-shellz(i_point+1,1);
        H11 = ysy;
        H2 = -shellz(i_point+1,1);
        H22 = shellz(i_point+1,2);
    end
    if ~ismember(i_point-1,long_point(:,1)) 
        if ~isnan(shelly(i_point+1,1)) && ~isnan(shellz(i_point+1,1)) 
            fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',H1,H11,z_height(c),H2,H22,z_height(c));
            fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
            fprintf(fid,'X %.3f S %.3f T %.3f\n',H1+(-shellz(i_point+1,1)-H2),H11,z_height(c));
            fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',H1+(-shellz(i_point+1,1)-H2),H11,z_height(c));
            fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',shelly(i_point+1,1)-shellz(i_point+1,1),shelly(i_point+1,2),z_height(c));
            X1 = shelly(i_point+1,1)-shellz(i_point+1,1);
            X2 = shelly(i_point+1,2);
            A1 = -shellz(i_point+1,1);
            A2 = shellz(i_point+1,2);
        elseif isnan(shelly(i_point+1,1)) && ~isnan(shellz(i_point+1,1)) 
            fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',H1,H11,z_height(c),H2,H22,z_height(c));
            fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
            fprintf(fid,'X %.3f S %.3f T %.3f\n',((H1)+(shelly(i_point+2,1)-shellz(i_point+1,1)))*0.5,(H11+shelly(i_point+2,2))*0.5,z_height(c));
            fprintf(fid,'X %.3f S %.3f T %.3f\n',shelly(i_point+2,1)-shellz(i_point+1,1),shelly(i_point+2,2),z_height(c));
            X1 = shelly(i_point+2,1)-shellz(i_point+1,1);
            X2 = shelly(i_point+2,2);
            A1 = -shellz(i_point+1,1);
            A2 = shellz(i_point+1,2);
        elseif ~isnan(shelly(i_point+1,1)) && isnan(shellz(i_point+1,1)) 
            fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',H1,H11,z_height(c),H2,H22,z_height(c));
            fprintf(fid,'A %.3f D %.3f E %.3f\n',-shellz(i_point+2,1),shellz(i_point+2,2),z_height(c));
            fprintf(fid,'X %.3f S %.3f T %.3f\n',H1+(-shellz(i_point+2,1)-H2),H11,z_height(c));
            fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',H1+(-shellz(i_point+2,1)-H2),H11,z_height(c));
            fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',shelly(i_point+1,1)-shellz(i_point+2,1),shelly(i_point+1,2),z_height(c));
            X1 = shelly(i_point+1,1)-shellz(i_point+2,1);
            X2 = shelly(i_point+1,2);
            A1 = -shellz(i_point+2,1);
            A2 = shellz(i_point+2,2);
        end
    elseif ismember(i_point-1,long_point(:,1)) 
        never = 1;
    end 
end
    end
end
fprintf(fid,'X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',X1,-20,50,A1,-20,50);
if u1 ~= u2
    flag = 1;
    for i_point = size(s,1) : size(t,1)
        pointx = t(i_point,1);
        pointy = t(i_point,2);
        if u1 < u2 
            if ~isnan(pointx)
                if flag
                    if i_point ~= size(s,1)
                        fprintf(fid,'A %.3f D %.3f E %.3f\n',A1,A2,z_height(c));
                        fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',A1,A2,z_height(c));
                        flag = 0;
                    end
                end
                if flag == 0
                    if i_point == size(s,1)
                        fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-pointx,pointy,z_height(c));
                        A1 = -pointx;
                    else
                        if ~isnan(pointx)
                        fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-pointx,pointy,z_height(c));
                        A1 = -pointx;
                        nnnn = 1;
                        end
                    end
                else
                    continue
                end
            else
                fprintf(fid,'A %.3f D %.3f E %.3f\n',-t(i_point+1,1),t(i_point+1,2),z_height(c));
                A1 = -t(i_point+1,1);
            end
        else 
            if ~isnan(pointx)
                if flag
                    if i_point ~= size(s,1)
                        fprintf(fid,'X %.3f S %.3f T %.3f\n',X1,X2,z_height(c));
                        fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',X1,X2,z_height(c));
                        flag = 0;
                    end
                end
                if flag == 0
                    if i_point == size(s,1)
                        fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',pointx+A1,pointy,z_height(c));
                        X1 = pointx + A1;
                    else
                        fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',pointx+A1,pointy,z_height(c));
                        X1 = pointx + A1;
                        nnnn = 0;
                    end
                else 
                    continue
                end
            else
                fprintf(fid,'X %.3f S %.3f T %.3f\n',t(i_point+1,1)+A1,t(i_point+1,2),z_height(c));
                X1 = t(i_point+1,1)+A1;
            end
        end
    end
    if nnnn == 1
        fprintf(fid,'A %.3f D %.3f E %.3f\n',A1,-20,50);
    elseif nnnn == 0
        fprintf(fid,'X %.3f S %.3f T %.3f\n',X1,-20,50);
    end
else
    time = 1;
end
elseif a == 0 
shelly = shell1;
shellz = shell2;
fprintf(fid,'X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',shelly(1,1)-shellz(1,1),shelly(1,2),z_height(c),-shellz(1,1),shellz(1,2),z_height(c));
u1 = length(shelly);
u2 = length(shellz);
if u1 > u2
    s = shellz;
    t = shelly;
else
    s = shelly;
    t = shellz;
end
flag = 1;
for i_point = 2 : size(s,1)
    zsx = shellz(i_point,1);
    zsy = shellz(i_point,2);
    ysx = shelly(i_point,1);
    ysy = shelly(i_point,2);
    if ~isnan(zsx) && ~isnan(ysx)
        if flag
            if i_point ~= 2
                fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                flag = 0;
            end
        end
        if flag == 0
            if i_point == 2
                fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                X1 = ysx-zsx;
                X2 = ysy;
                A1 = -zsx;
                A2 = zsy;
            else
                fprintf(fid,'G1 G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-zsx,ysy,z_height(c),-zsx,zsy,z_height(c));
                X1 = ysx-zsx;
                X2 = ysy;
                A1 = -zsx;
                A2 = zsy;
                never = -1;
            end
        else
            continue
        end
    elseif ~isnan(zsx) && isnan(ysx)
        fprintf(fid,'G2 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',shelly(i_point+1,1)-zsx,shelly(i_point+1,2),z_height(c),-zsx,zsy,z_height(c));
        X1 = shelly(i_point+1,1)-zsx;
        X2 = shelly(i_point+1,2);
        A1 = -zsx;
        A2 = zsy;
    elseif isnan(zsx) && ~isnan(ysx)
        fprintf(fid,'G1 X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',ysx-shellz(i_point+1,1),ysy,z_height(c),-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
        X1 = ysx-shellz(i_point+1,1);
        X2 = ysy;
        A1 = -shellz(i_point+1,1);
        A2 = shellz(i_point+1,2);
    else
        fprintf(fid,'X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',shelly(i_point+1,1)-shellz(i_point+1,1),shelly(i_point+1,2),z_height(c),-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
    end
end
fprintf(fid,'X %.3f S %.3f T %.3f A %.3f D %.3f E %.3f\n',X1,-20,50,A1,-20,50);
if u1 ~= u2
    flag = 1;
    for i_point = size(s,1) : size(t,1)
        pointx = t(i_point,1);
        pointy = t(i_point,2);
        if u1 < u2 
            if ~isnan(pointx)
                if flag
                    if i_point ~= size(s,1)
                        fprintf(fid,'A %.3f D %.3f E %.3f\n',A1,A2,z_height(c));
                        fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',A1,A2,z_height(c));
                        flag = 0;
                    end
                end
                if flag == 0
                    if i_point == size(s,1)
                        fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-pointx,pointy,z_height(c));
                        A1 = -pointx;
                    else
                        if ~isnan(pointx)
                        fprintf(fid,'G2 A %.3f D %.3f E %.3f\n',-pointx,pointy,z_height(c));
                        A1 = -pointx;
                        nnnn = 1;
                        end
                    end
                else
                    continue
                end
            else 
                fprintf(fid,'A %.3f D %.3f E %.3f\n',-t(i_point+1,1),t(i_point+1,2),z_height(c));
                A1 = -t(i_point+1,1);
            end
        else 
            if ~isnan(pointx)
                if flag
                    if i_point ~= size(s,1)
                        fprintf(fid,'X %.3f S %.3f T %.3f\n',X1,X2,z_height(c));
                        fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',X1,X2,z_height(c));
                        flag = 0;
                    end
                end
                if flag == 0
                    if i_point == size(s,1)
                        fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',pointx+A1,pointy,z_height(c));
                        X1 = pointx + A1;
                    else
                        fprintf(fid,'G1 X %.3f S %.3f T %.3f\n',pointx+A1,pointy,z_height(c));
                        X1 = pointx + A1;
                        nnnn = 0;
                    end
                else 
                    continue
                end
            else
                fprintf(fid,'X %.3f S %.3f T %.3f\n',t(i_point+1,1)+A1,t(i_point+1,2),z_height(c));
                X1 = t(i_point+1,1)+A1;
            end
        end
    end
    if nnnn == 1
        fprintf(fid,'A %.3f D %.3f E %.3f\n',A1,-20,50);
    elseif nnnn == 0
        fprintf(fid,'X %.3f S %.3f T %.3f\n',X1,-20,50);
    end
else
    time = 1;
end
end
end