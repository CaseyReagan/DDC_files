clear all
close all
clc

% fd=1;
% fs_org=8;
% L=256;
% ai=sign(randn(1,L));
% aq=sign(randn(1,L));
% U=fs_org/fd;
% a=ai+j*aq;
% b=upsample(a,U);
% fb=rcosfir(0.35,3,U,1,'sqrt');
% 
% c=filter(fb,1,b);
% 
% d=resample(c,1866,80);
% 
% e=d.*exp(j*2*pi*46.666/186.6*(1:length(d)));
% 
% f=real(e);

% load datasource.mat;

%load sig_source.dat

fid=fopen('C:\Users\Administrator\Desktop\Clk110_AD1_data.bin','r');
f = fread(fid,553648128,'int16');
% f = fread(fid,38000000,'int16');
fclose(fid);
%%
fs=110.0;
fd=1.333;
fc=10.0;
fs_new=8.0;

f = double(f);
f = f';
% f = f*10000;
% f = f*10;

gii=cos(2*pi*fc/fs*(1:length(f)));
gi=f.*cos(2*pi*fc/fs*(1:length(f)));

gq=f.*sin(2*pi*fc/fs*(1:length(f)));

fb1=fir1(128,4*fd/(fs/2));

hi=filter(fb1,1,gi);
hq=filter(fb1,1,gq);

%%
% y1 = 0;
% y2 = 0;
% y3 = 0;
% 
% x1 = 0;
% x2 = 0;
% x3 = 0;
% 
% for a=1:length(hi)
%     y1 = hi(a) + y1;
%     y2 = y1 + y2;
%     y3 = y2 + y3;
%     y4(a) = y3;
%     
%     x1 = hq(a) + x1;
%     x2 = x1 + x2;
%     x3 = x2 + x3;
%     x4(a) = x3;
% end
% 
% y5 = y4(1:8:end);
% x5 = x4(1:8:end);
% 
% x6 = 0;
% x7 = 0;
% x8 = 0;
% 
% y6 = 0;
% y7 = 0;
% y8 = 0;
% 
% ym1 = 0;
% ym2 = 0;
% ym3 = 0;
% 
% xm1 = 0;
% xm2 = 0;
% xm3 = 0;
% 
% for c=1:length(y5)
% %    for cs=1:comb_stage
%         ym1 = y5(c) - y6;
%         y6 = y5(c);
%         ym2 = ym1 - y7;
%         y7 = ym1;
%         ym3 = ym2 - y8;
%         y8 = ym2;
%         pi(c) = ym3;
%         
%         xm1 = x5(c) - x6;
%         x6 = x5(c);
%         xm2 = xm1 - x7;
%         x7 = xm1;
%         xm3 = xm2 - x8;
%         x8 = xm2;
%         pq(c) = xm3;
% end

% % my_cic_di=pi;
% % my_cic_dq=pq;
% % 
% pi = pi / 512;
% pq = pq / 512;
%%
y1 = 0;
y2 = 0;
y3 = 0;

x1 = 0;
x2 = 0;
x3 = 0;

for a=1:length(hi)
    y1 = hi(a) + y1;
    y2 = y1 + y2;
    y4(a) = y2;
    
    x1 = hq(a) + x1;
    x2 = x1 + x2;
    x4(a) = x2;
end

y5 = y4(1:8:end);
x5 = x4(1:8:end);

x6 = 0;
x7 = 0;
x8 = 0;

y6 = 0;
y7 = 0;
y8 = 0;

ym1 = 0;
ym2 = 0;
ym3 = 0;

xm1 = 0;
xm2 = 0;
xm3 = 0;

for c=1:length(y5)
%    for cs=1:comb_stage
        ym1 = y5(c) - y6;
        y6 = y5(c);
        ym2 = ym1 - y7;
        y7 = ym1;
        pi(c) = ym2;
        
        xm1 = x5(c) - x6;
        x6 = x5(c);
        xm2 = xm1 - x7;
        x7 = xm1;
        pq(c) = xm2;
end

% my_cic_di=pi;
% my_cic_dq=pq;
% 
% pi = pi / 512;
% pq = pq / 512;
%%
% hm = mfilt.cicdecim(8,1,3);
% 
% pi = filter(hm,hi);
% pq = filter(hm,hq);
% 
% pi=double(pi);
% pq=double(pq);



%%
td1 = 8/110.00;
td2 = 1/(fd*fs_new);

t1 = 1:td1:length(pi)*td1;
t2 = 1:td2:length(pi)*td1;

tl = td1;
d_t = t2(2) - t1(2);    % td2 - td1

c_1 = 0;
c0 = 0;
c1 = 0;
c2 = 0;

for d=1:length(t2)
    
%     if d == 1
%         u = 0.085;
%         j = 1;
%      end
    if d == 1
       u = 0.038;
       j = 1;
    end
    
    if d > 1
        u = u + d_t;
        if u > tl
            u = u - tl;
            j = j + 1;
        end
    end
%     yy(d) = u / tl * (pi(j+1) - pi(j)) + pi(j);
%     xx(d) = u / tl * (pq(j+1) - pq(j)) + pq(j);
    c_1 = 1/6 * u^3 - 1/6 * u;
    c0 = -1/2 * u^3 + 1/2 * u^2 + u;
    c1 = 1/2 * u^3 - u^2 - 1/2 * u + 1;
    c2 = -1/6 * u^3 + 1/2 * u^2 -1/3 * u;
    yy(d) = (c_1 * pi(j+3) + c0 * pi(j+2) + c1 * pi(j+1) + c2 * pi(j));
    xx(d) = (c_1 * pq(j+3) + c0 * pq(j+2) + c1 * pq(j+1) + c2 * pq(j));
    j = j + 1;
end

% for e=1:length(t2)
%     
%     if e == 1
%         u = 0;
%         j = 1;
%     else if e == 2
%             u = d_t;
%         end
%     end
%     
%     if e > 2
%         u = u + d_t;
%         if u > tl
%             u = u - tl;
%             j = j + 1;
%         end
%     end
%     xx(e) = u / tl * (pq(j+1) - pq(j)) + pq(j);
%     j = j + 1;
% end

my_intp_di=yy;
my_intp_dq=xx;

fb=rcosfir(0.35,3,fs_new/1,1,'sqrt');

qi=filter(fb,1,yy);
qq=filter(fb,1,xx);

my_fir_di=qi;
my_fir_dq=qq;

% zz=qi+sqrt(-1)*qq;

% my_zz = zz;
% % zz=pi+sqrt(-1)*pq;
% scatterplot(zz(4:8:end));
% scatterplot(zz(5:8:end));
% scatterplot(zz(6:8:end));

