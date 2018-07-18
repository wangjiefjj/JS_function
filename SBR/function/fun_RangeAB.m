function [ R ] = fun_RangeAB( Aw,Aj,Bw,Bj )
%http://blog.sina.com.cn/s/blog_658a93570101hynw.html
%%���������AB�����������룬Aw��γ�ȣ�Aj: ���� deg
Aw = Aw/180*pi;
Aj = Aj/180*pi;
Bw = Bw/180*pi;
Bj = Bj/180*pi;
Re = 6373e3;
c = acos(cos(pi/2-Bw)*cos(pi/2-Aw)+sin(pi/2-Bw)*sin(pi/2-Aw)*cos(Bj-Aj));
R = c*Re;
end

