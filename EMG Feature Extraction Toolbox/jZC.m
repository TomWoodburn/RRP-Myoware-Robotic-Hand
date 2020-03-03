%-------------------------------------------------------------------------%
%  Electromyography (EMG) Feature Extraction source codes demo version    %
%                                                                         %
%  Programmer: Jingwei Too                                                %
%                                                                         %
%  E-Mail: jamesjames868@gmail.com                                        %
%-------------------------------------------------------------------------%

%X = importdata('20_02_protocol/20_02_Protocol_light_tool_1.csv');

function ZC=jZC(X,thres)
N=length(X); ZC=0;
for i=1:N-1
  if ((X(i) > 0 && X(i+1) < 0) || (X(i) < 0 && X(i+1) > 0)) ...
      && (abs(X(i)-X(i+1)) >= thres)
    ZC=ZC+1;
  end
end
end

