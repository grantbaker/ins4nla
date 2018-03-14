a = dir('u0*.txt');
b = dir('v0*.txt');

for i = 1:length(a)
    u = load(a(i).name);   
    v = load(b(i).name);
    mesh(u.^2+v.^2)
    %axis([0 100 0 100 0 1e-6]) 
    drawnow
end
    