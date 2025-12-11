clear all
% load scaffold.mat
gfilename = 'Scaffold_without_MN-CAPP.txt';
fid = fopen(gfilename,'w');
%%
fprintf(fid,'X %.3f Y %.3f Z %.3f S %.3f T %.3f A %.3f B %.3f C %.3f D %.3f E %.3f\n',0,-20,50,20,50,0,-20,50,20,50);
z_height = 0.42 : 0.42 : 0.42 * 10;
for c = 1 : 10
    fprintf(fid,'Z %.3f C %.3f\n',z_height(c),z_height(c));
    shell1 = scaffold{c,1};
    shell2 = scaffold_small{c,1};
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
for i_point = 1 : size(s,1)
    zsx = shellz(i_point,1);
    zsy = shellz(i_point,2);
    ysx = shelly(i_point,1);
    ysy = shelly(i_point,2);
    if ~isnan(zsx) && ~isnan(ysx)
        if flag
            if i_point ~= 1
                fprintf(fid,'G1 G2 X %.3f Y %.3f Z %.3f A %.3f B %.3f C %.3f\n',shelly(1,1)-shellz(1,1),shelly(1,2),z_height(c),-shellz(1,1),shellz(1,2),z_height(c));
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
        fprintf(fid,'X %.3f Y %.3f Z %.3f A %.3f B %.3f C%.3f\n',shelly(i_point+1,1)-shellz(i_point+1,1),shelly(i_point+1,2),z_height(c),-shellz(i_point+1,1),shellz(i_point+1,2),z_height(c));
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