function make_intersections(s,ct,p)

user_input;
[nz,ny,nx]=size(s);

write=true;
% east
[k_east,r_east] = gamma_intersections(s,ct,p,-ny);
if ~zonally_periodic
    k_east(:,:,end)=nan;
    r_east(:,:,end)=nan;
end

% west
[k_west,r_west] = gamma_intersections(s,ct,p,ny);
if ~zonally_periodic
    k_west(:,:,1)=nan;
    r_west(:,:,1)=nan;
end

% north
[k_north,r_north] = gamma_intersections(s,ct,p,-1);
k_north(:,end,:)=nan;
r_north(:,end,:)=nan;

% south
[k_south,r_south] = gamma_intersections(s,ct,p,1);
k_south(:,1,:)=nan;
r_south(:,1,:)=nan;

if write   
    vars = {'k_east', 'r_east','k_west', 'r_west',...
        'k_north', 'r_north','k_south', 'r_south'};
    save('./data/intersections.mat', vars{:});  
end
    


