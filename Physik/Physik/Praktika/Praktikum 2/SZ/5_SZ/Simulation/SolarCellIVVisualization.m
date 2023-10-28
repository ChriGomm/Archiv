% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
% SolarCellIVVisualization HELP (This is short description)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ................
% ABOUT THE AUTHOR
% ................
%
% Name          : Zainul Aabdin
% Email         : zainul.aabdin@gmail.com
% About Me      : Research scholar (Ph.D.) in the group of 'Prof Oliver Eibl' at
%                 Eberhard Karls (University of Tuebingen), Tuebingen, Germany.
%                 You can know more about me by visiting my web page.
%
% ..................
% ABOUT THE PROGRAMM
% ..................
%
% Name          : SolarCellIVVisualization
% Version       : version 2.0 (February, 2014)
% Description   : This will simulate an I-V curve for solar cell
% 
% .................
% SPECIAL THANKS
% .................
% 
% Prof. Dr. Oliver Eibl         : I am very thankful to Prof Eibl, for
%                                 teaching me the first lesson of matlab and
%                                 motivating me to use matlab.
% Mr. Michael Dürrschnabel      : I am very thankful to Michael Durrschnabel for
%                                 his valuable suggestions.
%
% .................
% COPY RIGHT NOTICE (Fatima Website & Software Corporation @ All rights reserved)
% .................
%
% This script is free for students. They can use it, edit it and distribute
% it without any permission, but please don't remove the credit. Users are
% most welcome for comments, suggestions, bug reporting or any kind of help
% from me.
%
% Short Circuit Densities, Reverese Saturatution Current Densities, Band gap Energies
% and Temperature coeefficients from 
% Priyanka Singh and N.M.Ravindra, Solar Energy Materials & Solar Cells 100 (2012) 36-45
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% THANK YOU VERY MUCH TO ALL WHO USES THIS PROGRAMM...........REGARDS-ZAINUL
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function varargout = SolarCellIVVisualization(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SolarCellIVVisualization_OpeningFcn, ...
                   'gui_OutputFcn',  @SolarCellIVVisualization_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
end

function SolarCellIVVisualization_OpeningFcn(hObject, eventdata, handles, varargin)
clc
SolarIllumChange(hObject, eventdata, handles);
CellTempChange(hObject, eventdata, handles);
SerResChange(hObject, eventdata, handles);
ParResChange(hObject, eventdata, handles);
PlotFunction(handles);
movegui(handles.figure1,'center');
handles.output = hObject;
guidata(hObject, handles);
end

function varargout = SolarCellIVVisualization_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
end

% Ideal Plot Function~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% =========================================================================
function PlotFunction(handles)
global Iideal Uideal
S0=1367;             % Solar constant [W/m^2]
A0=1;                % Reference Area [cm^2]
k = 1.3806488e-23;   % Boltzman constant [J/K]
q = 1.6022e-19;	     % electronic change in [C]
Tnom = 298;          % Nominal temperature [K]

Eg0 = str2double(get(handles.EBandGapIdeal,'String'));              % Band gap of Si [eV] at 0 K
alpha_Eg=str2double(get(handles.EBandGapTempCoeffAlpha,'String'));  % band gap temperature coefficient [eV K^-1]
beta_Eg=str2double(get(handles.EBandGapTempCoeffBeta,'String'));    % band gap temperature coefficient [K]

A = str2double(get(handles.EExposedArea,'String'));                 % exposed area [cm^-2]
n = str2double(get(handles.EIdealiFacIdeal,'String'));              % diode ideality factor [unitless]
Isc = str2double(get(handles.EShortCircCurrIdeal,'String'))*10^-3;  % short circuit current density at 25 °C [A cm^-2]
Ki = str2double(get(handles.ETempCoeffIdeal,'String'))*10^-3;	    % short circuit current temperature coefficient [A cm^-2 K^-1]
s = str2double(get(handles.ERevSatTempExp,'String'));               % temperature exponent for diode reverse saturation current
I0nom = str2double(get(handles.ERevSatCurrIdeal,'String'))*10^-12;     % diode reverse saturation current density at 25 °C [A cm^-2]


S1 = str2double(get(handles.ESolarIllumIdeal,'String'));            % solar radiation in [W/m2]
T1 = str2double(get(handles.ECellTempIdeal,'String'));              % cell temperature [K]
Rs1 = str2double(get(handles.ESerResIdeal,'String'));               % cell series resitance [Ohm cm^2]
Rp1 = str2double(get(handles.EParResIdeal,'String'));              % cell shunt resitance or parallel resitance [Ohm cm^2]

S2 = str2double(get(handles.LSolarIllumCurrent,'String'));          % solar radiation in [W/m2]
T2 = str2double(get(handles.LCellTempCurrent,'String'));            % cell temperature [K]
Rs2 = str2double(get(handles.LSerResCurrent,'String'));             % cell series resitance [Ohm cm^2]
Rp2 = str2double(get(handles.LParResCurrent,'String'));            % cell shunt resitance or parallel resitance [Ohm cm^2]

U_D = 0:0.001:10;                                                   % Diode bias voltage [V]
cla(handles.axes1);

beta = (k*T1)/q;                                                    % diode thermal voltage [eV]
Eg=Eg0-alpha_Eg*T1^2/(beta_Eg+T1);                                  % band gap at current temperature [eV]
Iph = (Isc+Ki*(T1-298))*(S1/S0)*A;                                  % photo current [A]
I0 = I0nom*(T1/Tnom).^s*exp((Eg/1/beta)*((T1/Tnom)-1))*A;           % diode reverse saturation current at temperature T [A]
I=Iph-U_D/(Rp1/A)-I0*(exp(U_D/n/beta)-1);                           % cell current [A]
U=U_D-(Rs1/A)*I;                                                    % cell voltage [V]

h1 = plot(U,I*1000,'Color','b','LineWidth',2);
hold on;

Itest=I.*(I>=0);
indexIMin = find(Itest == min(Itest),1,'first');
VocIdeal = U(indexIMin);
set(handles.VocIdeal,'String', strcat(num2str(VocIdeal,'%.3f'),','));
set(handles.VocModified,'String', '---,');

Vtest=U.*(U>=0);
indexVMin = find(Vtest == min(Vtest),1,'last');
IscIdeal = I(indexVMin);
set(handles.IscIdeal,'String', strcat(num2str(IscIdeal*1000,'%.1f'),','));
set(handles.IscModified,'String', '---,');

step=ceil((indexIMin-indexVMin)/300);
indexstart=indexVMin-5*step;
indexstart=indexstart*(indexstart>=1)+indexVMin*(indexstart<1);
indexend=indexIMin+5*step;
indexend=indexend*(indexend<=numel(I))+indexIMin*(indexend>numel(I));
Iideal=I(indexstart:step:indexend)*1000;
Uideal=U(indexstart:step:indexend);

P = Vtest.*Itest;
indexMaxPow = find(P == max(P),1,'first');
PMPPIdeal = P(indexMaxPow);
set(handles.PMPPIdeal,'String', strcat(num2str(PMPPIdeal*1000,'%.1f'),','));
set(handles.PMPPModified,'String', '---,');

h3 = plot(U(indexMaxPow),I(indexMaxPow)*1000,'o','Color','k','Markersize',8,'MarkerFaceColor','b');

FFIdeal = (PMPPIdeal/(IscIdeal*VocIdeal))*100;
set(handles.FFIdeal,'String', strcat(num2str(FFIdeal,'%.1f'),','));
set(handles.FFModified,'String', '---,');

EtaIdeal=PMPPIdeal/S0/(A*1e-4)*100;
set(handles.EtaIdeal,'String', strcat(num2str(EtaIdeal,'%.1f'),','));
set(handles.EtaModified,'String', '---,');

indI_oc=find(I<=0,1,'first');
indI_p01Isc=find(I>=0.1*IscIdeal,1,'last');
indI_m01Isc=find(I<=-0.1*IscIdeal,1,'first');
RsIdeal=-(U(indI_m01Isc)-U(indI_p01Isc))/(I(indI_m01Isc)-I(indI_p01Isc))*A;
set(handles.RsIdeal,'String', strcat(num2str(RsIdeal,'%.1f'),','));
set(handles.RsModified,'String', '---,');

indV_sc=find(U<=0,1,'last');
indV_p01Voc=find(U>=0.1*VocIdeal,1,'first');
indV_m01Voc=find(U<=0.1*VocIdeal,1,'last');
RpIdeal=-(U(indV_p01Voc)-U(indV_m01Voc))/(I(indV_p01Voc)-I(indV_m01Voc))*A;
set(handles.RpIdeal,'String', strcat(num2str(RpIdeal,'%.1f'),''));
set(handles.RpModified,'String', '---');

leg = legend([h1,h3],'I-V Ideal','MPP Ideal');
legtxt = findobj(leg,'type','text');
%set(legtxt(1),'color','b');
%set(legtxt(2),'color','b');

grid on;

set(gca,'Xcolor',[0.5 0.5 0.5]);
set(gca,'Ycolor',[0.5 0.5 0.5]);
set(gca,'Color',[0.94 0.94 0.94]);

IscModified=0;
VocModified=0;
if (S1 ~= S2 || T1 ~= T2 || Rs1 ~= Rs2 || Rp1 ~= Rp2)
    beta = (k*T2)/q;                                                 % diode thermal voltage [eV]
    Eg=Eg0-alpha_Eg*T2^2/(beta_Eg+T2);                               % band gap at current temperature [eV]
    Iph = (Isc+Ki*(T2-298))*(S2/S0)*A;                               % photo current [A]
    I0 = I0nom*(T2/Tnom).^s*exp((Eg/1/beta)*((T2/Tnom)-1))*A;          % diode reverse saturation current at temperature T [A]
    I=U_D/(Rp2/A)+I0*(exp(U_D/n/beta)-1)-Iph;                      % cell current [A]
    U=U_D+(Rs2/A)*I;                                                 % cell voltage [V]
    I=-I;
    
    h2 = plot(U,I*1000,'Color','r','LineWidth',2);
    
    Itest=I.*(I>=0);
    indexIMin = find(Itest == min(Itest),1,'first');
    VocModified = U(indexIMin);
    set(handles.VocModified,'String', strcat(num2str(VocModified,'%.3f'),','));
    
    Vtest=U.*(U>=0);
    indexVMin = find(Vtest == min(Vtest),1,'last');
    IscModified = I(indexVMin);
    set(handles.IscModified,'String', strcat(num2str(IscModified*1000,'%.1f'),','));
    
    P = Vtest.*Itest;
    indexMaxPow = find(P == max(P),1,'first');
    PMPPModified = P(indexMaxPow);
    set(handles.PMPPModified,'String', strcat(num2str(PMPPModified*1000,'%.1f'),','));
    
    h4 = plot(U(indexMaxPow),I(indexMaxPow)*1000,'o','Color','k','Markersize',8,'MarkerFaceColor','r');
    
    FFModified = (PMPPModified/(IscModified*VocModified))*100;
    set(handles.FFModified,'String', strcat(num2str(FFModified,'%.1f'),','));
    
    EtaModified=PMPPModified/S0/(A*1e-4)*100;
    set(handles.EtaModified,'String', strcat(num2str(EtaModified,'%.1f'),','));
    
    indI_oc=find(I<=0,1,'first');
    indI_p01Isc=find(I>=0.1*IscModified,1,'last');
    indI_m01Isc=find(I<=-0.1*IscModified,1,'first');
    RsModified=-(U(indI_m01Isc)-U(indI_p01Isc))/(I(indI_m01Isc)-I(indI_p01Isc))*A;
    set(handles.RsModified,'String', strcat(num2str(RsModified,'%.1f'),','));
    
    indV_sc=find(U<=0,1,'last');
    indV_p01Voc=find(U>=0.1*VocModified,1,'first');
    indV_m01Voc=find(U<=0.1*VocModified,1,'last');
    RpModified=-(U(indV_p01Voc)-U(indV_m01Voc))/(I(indV_p01Voc)-I(indV_m01Voc))*A;
    set(handles.RpModified,'String', strcat(num2str(RpModified,'%.1f'),''));
        
    leg = legend([h1,h3,h2,h4],'I-V Ideal','MPP Ideal','I-V Modified','MPP Modified');
    legtxt = findobj(leg,'type','text');
    set(legtxt(1),'color','r');
    set(legtxt(2),'color','r');
    set(legtxt(3),'color','b');
    set(legtxt(4),'color','b');
end

x_min=0;
x_min_lin=x_min;
x_max=max(VocIdeal,VocModified);
x_max_lin=ceil(x_max);
x_tick_lin=[linspace(x_min_lin,x_max_lin,11)];

y_min=0;
y_min_lin=y_min;
y_max=max(IscIdeal,IscModified)*1000;
y_max_log=ceil(log10(y_max));
y_max_lin=10^(y_max_log);
y_tick_lin=[linspace(y_min_lin,y_max_lin,11)];

x_tick=x_tick_lin;
y_tick=y_tick_lin;
set(gca,'Units','normalized','Position',[0.15,0.1,0.8,0.85]);
set(gca,'Xtick',x_tick,'Ytick',y_tick,'Xlim',[min(x_tick) max(x_tick)],'Ylim',[min(y_tick) max(y_tick)]);
set(gca,'Xcolor','k','YColor','k','FontWeight','bold');
xlabel('Voltage (V)','Color','k','FontWeight','bold')
ylabel('Current (mA)','Color','k','FontWeight','bold')
set(leg,'color',[0.94 0.94 0.94]);
end

% Solar Illumination Functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% =========================================================================
function ESolarIllumIdeal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function ESolarIllumIdeal_Callback(hObject, eventdata, handles)
SolarIllumChange(hObject, eventdata, handles);
PlotFunction(handles);
end

function SSolarIllumChange_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end
function SSolarIllumChange_Callback(hObject, eventdata, handles)
SV = get(hObject,'Value');
set(handles.LSolarIllumCurrent,'String', num2str(SV,'%.1f'));
set(handles.LSolarIllumChange,'String', strcat(num2str((100*get(handles.SSolarIllumChange,'Value'))/str2double(get(handles.ESolarIllumIdeal,'String')),'%.1f'),' %'));
PlotFunction(handles);
end

function SolarIllumChange(hObject, eventdata, handles)
SV = str2double(get(handles.ESolarIllumIdeal,'String'));
set(handles.SSolarIllumChange,'max',2*SV,'min',0,'Value',SV)
set(handles.LSolarIllumCurrent,'String', num2str(SV,'%.1f'));
set(handles.LSolarIllumChange,'String', strcat(num2str((100*get(handles.SSolarIllumChange,'Value'))/SV,'%.1f'),' %'));
end


% Solar Cell Temperature Functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% =========================================================================
function ECellTempIdeal_Callback(hObject, eventdata, handles)
CellTempChange(hObject, eventdata, handles);
PlotFunction(handles);
end
function ECellTempIdeal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function SCellTempChange_Callback(hObject, eventdata, handles)
SV = get(hObject,'Value');
set(handles.LCellTempCurrent,'String', num2str(SV,'%.1f'));
set(handles.LCellTempChange,'String', strcat(num2str((100*get(handles.SCellTempChange,'Value'))/str2double(get(handles.ECellTempIdeal,'String')),'%.1f'),' %'));
PlotFunction(handles);
end
function SCellTempChange_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

function CellTempChange(hObject, eventdata, handles)
SV = str2double(get(handles.ECellTempIdeal,'String'));
set(handles.SCellTempChange,'max',2*SV,'min',0,'Value',SV)
set(handles.LCellTempCurrent,'String', num2str(SV,'%.1f'));
set(handles.LCellTempChange,'String', strcat(num2str((100*get(handles.SCellTempChange,'Value'))/SV,'%.1f'),' %'));
end


% Solar Cell Series Resistance Functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% =========================================================================
function ESerResIdeal_Callback(hObject, eventdata, handles)
SerResChange(hObject, eventdata, handles);
PlotFunction(handles);
end
function ESerResIdeal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function SSerResChange_Callback(hObject, eventdata, handles)
SV = get(hObject,'Value');
set(handles.LSerResCurrent,'String', num2str(SV,'%.1f'));
set(handles.LSerResChange,'String', strcat(num2str((100*get(handles.SSerResChange,'Value'))/str2double(get(handles.ESerResIdeal,'String')),'%.1f'),' %'));
PlotFunction(handles);
end
function SSerResChange_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

function SerResChange(hObject, eventdata, handles)
SV = str2double(get(handles.ESerResIdeal,'String'));
set(handles.SSerResChange,'max',50*SV,'min',0,'Value',SV)
set(handles.LSerResCurrent,'String', num2str(SV,'%.1f'));
set(handles.LSerResChange,'String', strcat(num2str((100*get(handles.SSerResChange,'Value'))/SV,'%.1f'),' %'));
end


% Solar Cell Parallel Resistance Functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% =========================================================================
function EParResIdeal_Callback(hObject, eventdata, handles)
ParResChange(hObject, eventdata, handles);
PlotFunction(handles);
end
function EParResIdeal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function SParResChange_Callback(hObject, eventdata, handles)
SV = get(hObject,'Value');
set(handles.LParResCurrent,'String', num2str(SV,'%.1f'));
set(handles.LParResChange,'String', strcat(num2str((100*get(handles.SParResChange,'Value'))/str2double(get(handles.EParResIdeal,'String')),'%.1f'),' %'));
PlotFunction(handles);
end
function SParResChange_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

function ParResChange(hObject, eventdata, handles)
SV = str2double(get(handles.EParResIdeal,'String'));
set(handles.SParResChange,'max',1.1*SV,'min',0,'Value',SV)
set(handles.LParResCurrent,'String', num2str(SV,'%.1f'));
set(handles.LParResChange,'String', strcat(num2str((100*get(handles.SParResChange,'Value'))/SV,'%.1f'),' %'));
end


% Solar Cell Short Circuit Current~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
% and Diode Reverse Saturartion Current properties~~~~~~~~~~~~~~~~~~~~~~~~~
% =========================================================================
function EShortCircCurrIdeal_Callback(hObject, eventdata, handles)
PlotFunction(handles);
end
function EShortCircCurrIdeal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function EIdealiFacIdeal_Callback(hObject, eventdata, handles)
PlotFunction(handles);
end
function EIdealiFacIdeal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ETempCoeffIdeal_Callback(hObject, eventdata, handles)
PlotFunction(handles);
end
function ETempCoeffIdeal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ERevSatCurrIdeal_Callback(hObject, eventdata, handles)
PlotFunction(handles);
end
function ERevSatCurrIdeal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function ERevSatTempExp_Callback(hObject, eventdata, handles)
PlotFunction(handles);
end
function ERevSatTempExp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function EExposedArea_Callback(hObject, eventdata, handles)
PlotFunction(handles);
end
function EExposedArea_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% Solar Cell Fixed Variables Functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% =========================================================================
function EBandGapTempCoeffBeta_Callback(hObject, eventdata, handles)
PlotFunction(handles);
end
function EBandGapTempCoeffBeta_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function EBandGapTempCoeffAlpha_Callback(hObject, eventdata, handles)
PlotFunction(handles);
end
function EBandGapTempCoeffAlpha_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function EBandGapIdeal_Callback(hObject, eventdata, handles)
PlotFunction(handles);
end
function EBandGapIdeal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% Save Screenshot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% =========================================================================
function Save_Screenshot_ClickedCallback(hObject, eventdata, handles)
Path='';
filename=strcat(Path,'Screenshot.jpg');
[Name,Path,Index] = uiputfile({'*.jpg';'*.*'},'Save Figure...',filename);
filename=strcat(Path,Name);
I = getframe(gcf);
imwrite(I.cdata, filename,'Quality',100);
PlotFunction(handles);
end



% Save I-V curve~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% =========================================================================
function Save_IV_curve_ClickedCallback(hObject, eventdata, handles)
Path='';
filename=strcat(Path,'I-V curve.jpg');
[Name,Path,Index] = uiputfile({'*.jpg';'*.*'},'Save Figure...',filename);
filename=strcat(Path,Name);
set(gcf,'Units','pixels');
rect=get(gcf,'Position');
I = getframe(gcf,[rect(3)*0.42,rect(4)*0.28,rect(3)*0.54,rect(4)*0.616]);
imwrite(I.cdata, filename,'Quality',100);
PlotFunction(handles);
end

% Save ideal I-V curve as Excel sheet~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% =========================================================================
function Save_ideal_IV_curve_as_Excel_sheet_ClickedCallback(hObject, eventdata, handles)
Path='';
filename=strcat(Path,'ideal I-V curve.xls');
[Name,Path,Index] = uiputfile({'*.xls';'*.*'},'Save ideal IV curve as Excel sheet...',filename);
filename=strcat(Path,Name);
global Iideal Uideal
lines=1;
IV_data{lines,1}='U (V)';
IV_data{lines,2}='I (mA)';
for lines=1:numel(Iideal);
IV_data{lines+1,1}=num2str(Uideal(lines),'%12.4f');
IV_data{lines+1,2}=num2str(Iideal(lines),'%12.4f');
end 
filename=strcat(Path,Name);
xlswrite(filename,IV_data);
PlotFunction(handles);
end

% Select Material~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% =========================================================================
function SelectMaterial_Callback(hObject, eventdata, handles)
%set(handles.SelectMaterial,'Value',0);
MaterialID=get(handles.SelectMaterial,'Value');
k = 1.3806488e-23;   % Boltzman constant [J/K]
q = 1.6022e-19;	     % electronic change in [C]
Tnom = 298;          % Nominal temperature [K]
beta = (k*Tnom)/q;   % diode thermal voltage [eV]
switch MaterialID
    case 1 
        set(handles.EBandGapIdeal,'String','1.1557');              % Band gap of Si [eV] at 0 K
        set(handles.EBandGapTempCoeffAlpha,'String','7.021e-04');   % band gap temperature coefficient [eV K^-1]
        set(handles.EBandGapTempCoeffBeta,'String','1108');        % band gap temperature coefficient [K]
        set(handles.ETempCoeffIdeal,'String','8.570e-03');          % short circuit current temperature coefficient [mA cm^-2 K^-1]
        set(handles.EShortCircCurrIdeal,'String','44.1');          % short circuit current density at 25 °C [mA cm^-2]
        set(handles.ERevSatTempExp,'String','0');                  % temperature exponent for diode reverse saturation current
        Eg0 = str2double(get(handles.EBandGapIdeal,'String'));              % Band gap [eV] at 0 K
        alpha_Eg=str2double(get(handles.EBandGapTempCoeffAlpha,'String'));  % band gap temperature coefficient [eV K^-1]
        beta_Eg=str2double(get(handles.EBandGapTempCoeffBeta,'String'));    % band gap temperature coefficient [K]
        Eg=Eg0-alpha_Eg*Tnom^2/(beta_Eg+Tnom);                              % band gap at nominal temperature [eV]
        m = str2double(get(handles.ERevSatTempExp,'String'));               % temperature exponent for diode reverse saturation current
        A=1.5e8;                                                   % Prefactor reverse saturation current density [mA cm^-2]
        I0nom = A*1e-3*Tnom^m*exp(-Eg/beta)*1e12;                     % diode reverse saturation current density at 25 °C [pA cm^-2]
        set(handles.ERevSatCurrIdeal,'String',num2str(I0nom,'%.3e')); % diode reverse saturation current density at 25 °C [pA cm^-2]
    case 2
        set(handles.EBandGapIdeal,'String','0.7412');              % Band gap of Ge [eV] at 0 K
        set(handles.EBandGapTempCoeffAlpha,'String','4.561e-04');   % band gap temperature coefficient [eV K^-1]
        set(handles.EBandGapTempCoeffBeta,'String','210');         % band gap temperature coefficient [K]
        set(handles.ETempCoeffIdeal,'String','0.106e-03');          % short circuit current temperature coefficient [mA cm^-2 K^-1]
        set(handles.EShortCircCurrIdeal,'String','61.0');          % short circuit current density at 25 °C [mA cm^-2]
        set(handles.ERevSatTempExp,'String','0');                  % temperature exponent for diode reverse saturation current
        Eg0 = str2double(get(handles.EBandGapIdeal,'String'));              % Band gap [eV] at 0 K
        alpha_Eg=str2double(get(handles.EBandGapTempCoeffAlpha,'String'));  % band gap temperature coefficient [eV K^-1]
        beta_Eg=str2double(get(handles.EBandGapTempCoeffBeta,'String'));    % band gap temperature coefficient [K]
        Eg=Eg0-alpha_Eg*Tnom^2/(beta_Eg+Tnom);                              % band gap at nominal temperature [eV]
        m = str2double(get(handles.ERevSatTempExp,'String'));               % temperature exponent for diode reverse saturation current
        A=1.5e8;                                                   % Prefactor reverse saturation current density [mA cm^-2]
        I0nom = A*1e-3*Tnom^m*exp(-Eg/beta)*1e12;                     % diode reverse saturation current density at 25 °C [pA cm^-2]
        set(handles.ERevSatCurrIdeal,'String',num2str(I0nom,'%.3e')); % diode reverse saturation current density at 25 °C [pA cm^-2]
    case 3
        set(handles.EBandGapIdeal,'String','1.5216');              % Band gap of GaAs [eV] at 0 K
        set(handles.EBandGapTempCoeffAlpha,'String','8.871e-04');   % band gap temperature coefficient [eV K^-1]
        set(handles.EBandGapTempCoeffBeta,'String','572');         % band gap temperature coefficient [K]
        set(handles.ETempCoeffIdeal,'String','0.196e-03');          % short circuit current temperature coefficient [mA cm^-2 K^-1]
        set(handles.EShortCircCurrIdeal,'String','31.6');          % short circuit current density at 25 °C [mA cm^-2]
        set(handles.ERevSatTempExp,'String','0');                  % temperature exponent for diode reverse saturation current
        Eg0 = str2double(get(handles.EBandGapIdeal,'String'));              % Band gap [eV] at 0 K
        alpha_Eg=str2double(get(handles.EBandGapTempCoeffAlpha,'String'));  % band gap temperature coefficient [eV K^-1]
        beta_Eg=str2double(get(handles.EBandGapTempCoeffBeta,'String'));    % band gap temperature coefficient [K]
        Eg=Eg0-alpha_Eg*Tnom^2/(beta_Eg+Tnom);                              % band gap at nominal temperature [eV]
        m = str2double(get(handles.ERevSatTempExp,'String'));               % temperature exponent for diode reverse saturation current
        A=1.5e8;                                                   % Prefactor reverse saturation current density [mA cm^-2]
        I0nom = A*1e-3*Tnom^m*exp(-Eg/beta)*1e12;                     % diode reverse saturation current density at 25 °C [pA cm^-2]
        set(handles.ERevSatCurrIdeal,'String',num2str(I0nom,'%.3e')); % diode reverse saturation current density at 25 °C [pA cm^-2]
    case 4
        set(handles.EBandGapIdeal,'String','1.4206');              % Band gap of InP [eV] at 0 K
        set(handles.EBandGapTempCoeffAlpha,'String','4.906e-04');   % band gap temperature coefficient [eV K^-1]
        set(handles.EBandGapTempCoeffBeta,'String','93');          % band gap temperature coefficient [K]
        set(handles.ETempCoeffIdeal,'String','9.430e-03');          % short circuit current temperature coefficient [mA cm^-2 K^-1]
        set(handles.EShortCircCurrIdeal,'String','34.7');          % short circuit current density at 25 °C [mA cm^-2]
        set(handles.ERevSatTempExp,'String','3');                  % temperature exponent for diode reverse saturation current
        Eg0 = str2double(get(handles.EBandGapIdeal,'String'));              % Band gap [eV] at 0 K
        alpha_Eg=str2double(get(handles.EBandGapTempCoeffAlpha,'String'));  % band gap temperature coefficient [eV K^-1]
        beta_Eg=str2double(get(handles.EBandGapTempCoeffBeta,'String'));    % band gap temperature coefficient [K]
        Eg=Eg0-alpha_Eg*Tnom^2/(beta_Eg+Tnom);                              % band gap at nominal temperature [eV]
        m = str2double(get(handles.ERevSatTempExp,'String'));               % temperature exponent for diode reverse saturation current
        C=17.9;                                                      % Prefactor reverse saturation current density [mA cm^-2 K-3]
        I0nom = C*1e-3*Tnom^m*exp(-Eg/beta)*1e12;                     % diode reverse saturation current density at 25 °C [pA cm^-2]
        set(handles.ERevSatCurrIdeal,'String',num2str(I0nom,'%.3e')); % diode reverse saturation current density at 25 °C [pA cm^-2]
    case 5
        set(handles.EBandGapIdeal,'String','1.6077');              % Band gap of CdTe [eV] at 0 K
        set(handles.EBandGapTempCoeffAlpha,'String','3.100e-04');   % band gap temperature coefficient [eV K^-1]
        set(handles.EBandGapTempCoeffBeta,'String','108');         % band gap temperature coefficient [K]
        set(handles.ETempCoeffIdeal,'String','0.121e-03');          % short circuit current temperature coefficient [mA cm^-2 K^-1]
        set(handles.EShortCircCurrIdeal,'String','28.2');          % short circuit current density at 25 °C [mA cm^-2]
        set(handles.ERevSatTempExp,'String','3');                  % temperature exponent for diode reverse saturation current
        Eg0 = str2double(get(handles.EBandGapIdeal,'String'));              % Band gap [eV] at 0 K
        alpha_Eg=str2double(get(handles.EBandGapTempCoeffAlpha,'String'));  % band gap temperature coefficient [eV K^-1]
        beta_Eg=str2double(get(handles.EBandGapTempCoeffBeta,'String'));    % band gap temperature coefficient [K]
        Eg=Eg0-alpha_Eg*Tnom^2/(beta_Eg+Tnom);                              % band gap at nominal temperature [eV]
        m = str2double(get(handles.ERevSatTempExp,'String'));               % temperature exponent for diode reverse saturation current
        C=17.9;                                                      % Prefactor reverse saturation current density [mA cm^-2 K-3]
        I0nom = C*1e-3*Tnom^m*exp(-Eg/beta)*1e12;                     % diode reverse saturation current density at 25 °C [pA cm^-2]
        set(handles.ERevSatCurrIdeal,'String',num2str(I0nom,'%.3e')); % diode reverse saturation current density at 25 °C [pA cm^-2]
    case 6
        set(handles.EBandGapIdeal,'String','2.583');               % Band gap of CdS [eV] at 0 K
        set(handles.EBandGapTempCoeffAlpha,'String','4.020e-04');   % band gap temperature coefficient [eV K^-1]
        set(handles.EBandGapTempCoeffBeta,'String','147');         % band gap temperature coefficient [K]
        set(handles.ETempCoeffIdeal,'String','6.670e-03');          % short circuit current temperature coefficient [mA cm^-2 K^-1]
        set(handles.EShortCircCurrIdeal,'String','7.5');           % short circuit current density at 25 °C [mA cm^-2]
        set(handles.ERevSatTempExp,'String','3');                  % temperature exponent for diode reverse saturation current
        Eg0 = str2double(get(handles.EBandGapIdeal,'String'));              % Band gap [eV] at 0 K
        alpha_Eg=str2double(get(handles.EBandGapTempCoeffAlpha,'String'));  % band gap temperature coefficient [eV K^-1]
        beta_Eg=str2double(get(handles.EBandGapTempCoeffBeta,'String'));    % band gap temperature coefficient [K]
        Eg=Eg0-alpha_Eg*Tnom^2/(beta_Eg+Tnom);                              % band gap at nominal temperature [eV]
        m = str2double(get(handles.ERevSatTempExp,'String'));               % temperature exponent for diode reverse saturation current
        C=17.9;                                                      % Prefactor reverse saturation current density [mA cm^-2 K-3]
        I0nom = C*1e-3*Tnom^m*exp(-Eg/beta)*1e12;                     % diode reverse saturation current density at 25 °C [pA cm^-2]
        set(handles.ERevSatCurrIdeal,'String',num2str(I0nom,'%.3e')); % diode reverse saturation current density at 25 °C [pA cm^-2]
end    
PlotFunction(handles);
end


% --- Executes during object creation, after setting all properties.
function SelectMaterial_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
