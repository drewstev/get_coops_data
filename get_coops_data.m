function data =  get_coops_data(site,varargin)
% GET_COOPS_DATA - retreive NOAA COOPs data using web services
%
%   DATA = GET_COOPS_DATA(SITE) - retrieves NOAA COOPS water level 
%       data for the station specified in SITE.  SITE should be a character
%       string (eg. '9447130' for Seattle) specifying a single station. 
%       The user will be prompted for the beginning and ending dates 
%       to download. The data are returned in the structure array, DATA. 
%
%   OPTIONAL INPUTS
%       Optional inputs are supplied as property pair or as an input
%       structure.
%
%       'start_time' - beginning of requested time-series (string)
%       'end_time'   -  end of requested time-series (string)
%       'parameter'  - type of data requested. The list of supported
%                      parameters includes: 
%                           'water_level' (default)
%                           'air_temperature'
%                           'water_temperature'
%                           'wind'
%                           'air_pressure'
%                           'conductivity'
%                           'humidity'
%                           'salinity'
%                           'hourly_height'
%                           'high_low'
%                           'daily_mean'
%                           'monthly_mean'
%                           'one_minute_water_level'
%                           'predictions'
%                           'datums'
%        'datum'     - Datum of requested output data. See the station 
%                      home page for a list of available datums at your 
%                      site.
%        'time_zone' - Time zone of requested output. Function supports
%                      'GMT' (default), 'LST', or 'LST_LDT'.
%        'units'     - Units of output. Either 'METRIC' (default) or 
%                      'ENGLISH'.
%        'interval'  - Not necessary for most requests. Set to 'hilo' 
%                      to retreive tide predictions at subordinate
%                      stations.
%                        
%
%   EXAMPLES
%       Example 1. User supplies the station id (required) and is prompted
%       for the start and end times.
%           site = '9447130';
%           data = get_coops_data(site)
%
%       Example 2. User requests wind data rather than water levels
%           site = '9447130';
%           data = get_coops_data(site,'parameter','wind')
%
%       Example 3. User supplies start and end time through property value
%       pairs. Get last week of data.
%           site = '9447130';
%           st = datestr(now-7);
%           et = datestr(now);
%           data = get_coops_data(site,'start_time',st,'end_time',et)
%
%       Example 4. User supplies optional inputs as a structure.
%           site = '9447130';
%           opt.start_time = datestr(now-7);
%           opt.end_time = datestr(now);
%           opt.datum = 'navd';
%           data = get_coops_data(site,opt);
%
%          %plot the results
%           figure,
%           plot(data.mtime,data.water_level)
%           ylabel(['Water level (m, ',opt.datum,')'])
%
%           datetick('x')
%
%       Example 5. Download the data for the last week for a group of sites
%       in the Puget Sound area
%           sites={'9443090';... %Neah Bay
%                  '9444090';... %Port Angeles
%                  '9444900';... %Port Townsend
%                  '9447130';... %Seattle
%                  '9446484'};   %Tacoma
%           opt.datum='MLLW';
%           opt.start_time = datestr(now-7);
%           opt.end_time = datestr(now);
%           data=cellfun(@(x)(get_coops_data(x,opt)),sites);
%
%           figure
%           hold on 
%           arrayfun(@(x)(plot(x.mtime,x.water_level)),data)
%           datetick('x')
%           ylabel(['Water level (m, ',opt.datum,')'])
%
%           stn_names=arrayfun(@(x)(x.site_name),data,'un',0);
%           legend(stn_names,'location','eastoutside')
%
%        Example 6. Plot predicted versus observed water levels for
%        Seattle during the Hanukkah Eve storm in 2006.
%           site = '9447130';
%           opt.datum='MSL';
%           opt.start_time = '2006-12-13';
%           opt.end_time = '2006-12-17';
%     
%           pred = get_coops_data(site,opt,'parameter','predictions');
%           obs  = get_coops_data(site,opt,'parameter','water_level');
%           wind = get_coops_data(site,opt,'parameter','wind');
%           
%           figure
%           subplot(211)
%           plot(pred.mtime,pred.wl_predictions)
%           hold on 
%           plot(obs.mtime,obs.water_level)     
%           datetick('x')
%           ylabel(['Water level (m, ',opt.datum,')'])
%           
%           %caluclate residual
%           wl_interp=interp1(obs.mtime,obs.water_level,pred.mtime);
%           wl_diff=wl_interp-pred.wl_predictions;
%           plot(pred.mtime,wl_diff)
%   
%           legend('Predicted','Observed','Residual')
%   
%           subplot(212)
%           plot(wind.mtime,wind.wind_speed)
%           datetick('x')
%           ylabel('Wind Speed (m/s)')
%
%        Example 7. Plot monthly mean sea level trends for Seattle and 
%        San Francisco between 1900 and today.
%           sites={'9414290';... % San Francisco
%                  '9447130'};   %Seattle
%           opt.start_time = '1900-1-1';
%           opt.end_time = datestr(now);
%           opt.parameter = 'monthly_mean';
%           opt.datum = 'MSL';
%           data=cellfun(@(x)(get_coops_data(x,opt)),sites);
%
%           figure
%           hold on 
%           arrayfun(@(x)(plot(x.mtime,x.MSL)),data)
%           datetick('x')
%           ylabel(['Water level (m, ',opt.datum,')'])
%
%           stn_names=arrayfun(@(x)(x.site_name),data,'un',0);
%           legend(stn_names,'location','northwest',...
%               'orientation','horiz')
%          
% More information on NOAA COOPS data, including a station listing,
% can be found <a href = "https://tidesandcurrents.noaa.gov/">here</a>.
% Information about access to NOAA COOPS data data products through web
% services is available <a href = "https://tidesandcurrents.noaa.gov/api/">here</a>.
%
% Requires Matlab 2016b or later 

% Andrew Stevens, astevens@usgs.gov
% version 1.02
% 12/11/2018


p=inputParser;
p.KeepUnmatched=1;

addRequired(p,'site',@ischar);
opts={'start_time', [],      {'char'},               {};...
    'end_time',   [],      {'char'},               {};...
    'parameter', 'water_level',      {'char'},     {};...
    'datum',     'MLLW',      {'char'},            {};...
    'units',     'metric',      {'char'},          {};...
    'time_zone', 'gmt',         {'char'},          {};...
    'interval',   [],           {'char'}        {}};

cellfun(@(x)(p.addParameter(x{1},x{2},...
    @(y)(validateattributes(y, x{3},x{4})))),num2cell(opts,2))
p.parse(site,varargin{:});
opt=p.Results;

params= {'water_level';...
    'air_temperature';...
    'water_temperature';...
    'wind';...
    'air_pressure';...
    'conductivity';...
    'humidity';...
    'salinity';...
    'hourly_height';...
    'high_low';...
    'daily_mean';...
    'monthly_mean';...
    'one_minute_water_level';...
    'predictions';...
    'datums'};

wlparams=params([1,9:14]);

validatestring(opt.parameter,....
    params,...
    'get_coops_data','parameter');
validatestring(opt.datum,...
    {'CRD';'IGLD';'LWD';'MHHW';'MHW';...
    'MTL';'MSL';'MLW';'MLLW';'NAVD';...
    'STND'},...
    'get_coops_data','datum');
validatestring(opt.units,....
    {'metric';'english'},...
    'get_coops_data','units');
validatestring(opt.time_zone,....
    {'gmt';'lst';'lst_ldt'},...
    'get_coops_data','time_zone');
if ~isempty(opt.interval)
    validatestring(opt.interval,....
        {'hilo';'h'},...
        'get_coops_data','interval');
end
%rename variable names to more user-friendly values
%also need to specify weather output is a string or number
switch opt.parameter
    case 'water_level'
        ds={'v'};
        dt={'n'};
        dval={'water_level'};
        meta=['Preliminary or verified water levels, ',...
            'depending on availability'];
        range=31;
    case 'wind'
        ds={'s';'d';'g'};
        dt={'n';'n';'n'};
        dval={'wind_speed';'wind_dir';'wind_gust'};
        meta=['Wind speed, direction, and gusts as measured ',...
            'at the station'];
        range=31;
    case 'air_temperature'
        ds={'v'};
        dt={'n'};
        dval={'air_temperature'};
        meta='Air temperature as measured at the station';
        range=31;
    case 'water_temperature'
        ds={'v'};
        dt={'n'};
        dval={'water_temperature'};
        meta='Water temperature as measured at the station';
        range=31;
    case 'air_pressure'
        ds={'v'};
        dt={'n'};
        dval={'air_pressure'};
        meta='Barometric pressure as measured at the station';
        range=31;
    case 'conductivity'
        ds={'v'};
        dt={'n'};
        dval={'conductivity'};
        meta='The water''s conductivity as measured at the station';
        range=31;
    case 'humidity'
        ds={'v'};
        dt={'n'};
        dval={'humidity'};
        meta='Relative humidity as measured at the station';
        range=31;
    case 'salinity'
        ds={'v'};
        dt={'n'};
        dval={'salinity'};
        meta='Salinity and specific gravity data for the station';
        range=31;
    case 'hourly_height'
        ds={'v'};
        dt={'n'};
        dval={'hourly_height'};
        meta='Verified hourly height water level data for the station';
        range=365;
    case 'high_low'
        ds={'v';'ty'};
        dt={'n';'f'};
        dval={'high_low_height';'high_low_type'};
        meta='Verified high/low water level data for the station';
        range=365;
    case 'daily_mean'
        ds={'v'};
        dt={'n'};
        dval={'daily_mean'};
        meta='Verified daily mean water level data for the station';
        range=365;
    case 'monthly_mean'
        ds={'highest';...
            'MHHW';...
            'MHW';...
            'MSL';...
            'MTL';...
            'MLW';...
            'MLLW';...
            'DTL';...
            'GT';...
            'MN';...
            'DHQ';...
            'DLQ';...
            'HWI';...
            'LWI';...
            'lowest'};
        dt=repmat({'n'},length(ds),1);
        dval=ds;
        meta='Verified monthly mean water level data for the station';
        range=365*10;
        
    case 'one_minute_water_level'
        ds={'v'};
        dt={'n'};
        dval={'one_minute_water_level'};
        meta='One minute water level data for the station';
        range=31;
        
    case 'predictions'
        ds={'v'};
        dt={'n'};
        dval={'wl_predictions'};
        meta='6 minute predictions water level data for the station';
        range=31;
    case 'datums'
        ds={'MHHW';...
            'MHW';...
            'DTL';...
            'MTL';...
            'MSL';...
            'MLW';...
            'MLLW';...
            'GT';...
            'MN';...
            'DHQ';...
            'DLQ';...
            'NAVD';...
            'LWI';...
            'HWI'};
        dt=repmat({'n'},length(ds),1);
        dval=ds;
        meta='datums data for the stations';
        range=365*10;
end

%deal with special case of high_low predictions
if strcmpi(opt.interval,'hilo')
        ds={'v';'type'};
        dt={'n';'f'};
        dval={'high_low_height';'high_low_type'};
        meta='High/low water level data for the station';
end

%start and end times
if isempty(opt.start_time)
    fprintf('Pick a start time.\n');
    dnstart=uigetdate(now);
else
    dnstart=datenum(opt.start_time);
end
if isempty(opt.end_time)
    fprintf('Pick an end time.\n');
    dnend=uigetdate(dnstart);
else
    dnend=datenum(opt.end_time);
end


if dnend-dnstart>range
    dnv=[dnstart:range:dnend,dnend];
    dnstart=dnv(1:end-1);
    dnend=dnv(2:end);
end

mtime=cell(length(dnstart),1);
vals=cell(length(dnstart),length(ds));

webopts=weboptions('Timeout',20);
for i=1:length(dnstart)
    
    url=['https://tidesandcurrents.noaa.gov/api/datagetter?',...
        'begin_date=',datestr(dnstart(i),'yyyymmdd HH:MM'),...
        '&end_date=',datestr(dnend(i),'yyyymmdd HH:MM'),...
        '&station=',site,...
        '&product=',opt.parameter,...
        '&datum=',opt.datum,...
        '&units=',opt.units,...
        '&time_zone=',opt.time_zone,...
        '&application=web_services&format=json'];
    
    if ~isempty(opt.interval)
        switch opt.interval
            case 'hilo'
                url=[url,'&interval=hilo']; %#ok
            case 'h'
                url=[url,'&interval=h']; %#ok
        end
    end
    
    rdata=jsondecode(webread(url,webopts));
    if isfield(rdata,'error')
        error(rdata.error)
    end
    if strcmpi(opt.parameter,'predictions')
        rdata.data=rdata.predictions;
    end
    
    %only allow one station per request
    if i==1
        if isfield(rdata,'metadata')
            data.site_name=rdata.metadata.name;
            data.site_id=rdata.metadata.id;
            data.site_lon=rdata.metadata.lon;
            data.site_lat=rdata.metadata.lat;
            data.time_zone=opt.time_zone;
            data.units=opt.units;
            data.description=meta;
            if ~isempty(intersect(opt.parameter,wlparams))
                data.datum=opt.datum;
            end
        else
           
            data.site_name=opt.site;
            data.time_zone=opt.time_zone;
            data.units=opt.units;
            data.description=meta;
            if ~isempty(intersect(opt.parameter,wlparams))
                data.datum=opt.datum;
            end
        end
    end
    
    if strcmpi(opt.parameter,'datums')
        for j=1:length(rdata.datums)
            data.(rdata.datums(j).n)=str2double(rdata.datums(j).v);
        end
%         data=rmfield(data,'time_zone');
        
    else
        
        for j=1:length(ds)
            if strcmpi(dt{j},'n')
                vals{i,j}=arrayfun(@(x)(str2double(x.(ds{j}))),...
                    rdata.data);
            else
                vals{i,j}=arrayfun(@(x)(x.(ds{j})),...
                    rdata.data,'un',0);
            end
        end
        
        %convert to datetimes
        if strcmpi(opt.parameter,'monthly_mean')
            yr=arrayfun(@(x)(str2double(x.year)),...
                rdata.data);
            mm=arrayfun(@(x)(str2double(x.month)),...
                rdata.data);
            mtime{i}=datenum(yr,mm,ones(length(yr),1));
            
        else
            dstr=arrayfun(@(x)(x.t),...
                rdata.data,'un',0);
            mtime{i}=datenum(dstr,'yyyy-mm-dd HH:MM');
        end
    end
end

if ~strcmpi(opt.parameter,'datums')
    mt=cell2mat(mtime);
    
    [data.mtime,midx]=unique(mt);
    for j=1:length(ds)
        if strcmpi(dt{j},'n')
            dm=cell2mat(vals(:,j));
        else
            dm=cat(1,vals{:,j});
        end
        data.(dval{j})=dm(midx,1);
    end
end
%-------------------------------------------------------------------------
function out = uigetdate(varargin)
% UIGETDATE  date selection dialog box
%    T = UIGETDATE(D) displays a dialog box in form of a calendar 
%    
%    UIGETDATE expects serial date number or standard MATLAB Date 
%    format (see DATESTR) as input data und returns serial date number 
%    for the selected date and time.
%
%    UIGETDATE by itself uses the current date and time as input data
%
% Example:
%         t = datestr( uigetdate('16-Aug-1974 03:00') )
% 
% See also datevec, datestr, datenum

%   version: v1.0
%   author:  Elmar Tarajan [MCommander@gmx.de]

if nargin == 0
   varargin{1} = now;
end% if

if ~ishandle(varargin{1})
   %
   datvec = datevec(varargin{1});
   %
   scr = get(0,'ScreenSize');
   h.units = 'pixels';
   h.parent = figure(h,'menubar','none', ...
            'numbertitle','off', ...
            'resize','off', ...
            'handlevisibility','on', ...
            'visible','off', ...            
            'WindowStyle','modal', ...
            'Tag','uigetdate', ...
            'position',[ (scr(3:4)- [197 199])/2 197 199 ]);
   %
   pos = [5 5 0 0];
   uicontrol(h,'style','edit','position',pos+[0 0 104 26])
   uicontrol('style','slider','units','pixels','position',pos+[3 2 100 20], ...
             'sliderstep',[.0005 .0005],'min',-10,'max',10,'value',0, ...
             'callback','uigetdate(gcbo,''time'')','UserData',0)
   %
   h.style           = 'edit';
   h.fontweight      = 'bold';
   h.foregroundcolor = [.2 .2 .2];
   uicontrol(h,'enable','inactive','position',pos+[ 17 2 73 20],'Tag','time', ...
               'String',sprintf('%02d:%02d',datvec(4:5)))
   %
   % textbanners
   tmp = [2 20 101 4 ; 2 2 101 3 ; 2 2 3 22 ; 17 2 2 22 ; 88 2 2 22 ; 101 2 2 22 ];
   for i=1:6 ; uicontrol(h,'style','text','position',pos+tmp(i,:)) ; end% for
   %
   uicontrol(h,'style','edit','position',pos+[105 0 84 26],'visible','on')   
   uicontrol(h,'style','pushbutton','position',pos+[108 2 78 21],'Tag','ok', ...
               'CData',repmat(repmat([0.3:0.01:1 1:-0.01:0.3],18,1),[1 1 3]), ...
               'string','ok','Callback','uigetdate(gcbo,''ok'')')
   %
   pos = [5 32 0 0];
   uicontrol(h,'style','edit','position',pos+[0 0 189 136],'enable','inactive','Tag','cday', ...
      'UserData',datvec(3))   
   h.style      = 'pushbutton';
   h.fontweight = 'normal';
   for i=95:-19:0
      for j=0:26:156
         uicontrol(h,'position',pos+[j+3 i+2 27 20],'Enable','off', ...
                     'foregroundcolor',[.2 .2 .2],'Tag','day', ...
                     'callback','uigetdate(gcbo,''day'')')
      end% for
   end% for
   %
   tmp = {'Mon' 'Tue' 'Wed' 'Thu' 'Fri' 'Sat' 'Sun'};
   for j=0:6
      uicontrol(h,'style','text','position',pos+[j*26+4 119 25 13],'string',tmp{j+1}, ...
                  'backgroundcolor',[0.4 0.4 0.4],'foregroundcolor',[.9 .9 .9])         
   end% for
   %
   pos = [5 169 0 0];
   uicontrol(h,'style','edit','position',pos+[0 0 189 26])
   h.style = 'slider';
   uicontrol(h,'position',pos+[3 2 100 20],'sliderstep',[0.00025 1], ...
               'min',-2000,'max',2000,'Value',datvec(2), ...
               'callback','uigetdate(gcbo,''months'')')
   uicontrol(h,'position',pos+[112 2 74 20],'sliderstep',[0.00025 1], ...
               'min',0,'max',4000,'value',datvec(1), ...
               'callback','uigetdate(gcbo,''year'')')
   %
   h.style           = 'edit';
   h.enable          = 'inactive';
   h.fontweight      = 'bold';
   h.foregroundcolor = [.2 .2 .2];
   tmp = {'Januar' 'Februar' 'March' 'April' 'May' 'Juni' 'Juli' ...
          'August' 'September' 'October' 'November' 'December'};
   uicontrol(h,'position',pos+[ 17 2 73 20],'Tag','months','String',tmp{datvec(2)},'Userdata',tmp)
   uicontrol(h,'position',pos+[126 2 47 20],'Tag','year','String',num2str(datvec(1)))
   %
   % textbanners
   h.style = 'text';
   tmp = [2 20 185 4 ; 2 2 185 3 ; 2 2 3 22 ; 17 2 2 22 ; 88 2 2 22 ; ...
      101 2 13 22 ; 126 2 2 22 ; 171 2 2 22 ; 184 2 3 22];
   for i=1:9
      uicontrol(h,'position',pos+tmp(i,:))
   end% for
   %
   set(h.parent,'visible','on')
   setday(varargin{1})
   %
   set(findobj(gcf,'string',num2str(datvec(3))),'CData',geticon)
   %
   uiwait
   try
      out = datenum([num2str( ...
               get(findobj(gcf,'Tag','cday'),'UserData')) '-' ...
               get(findobj(gcf,'Tag','months'),'String') '-' ...
               get(findobj(gcf,'Tag','year'),'String') ' ' ...
               get(findobj(gcf,'Tag','time'),'String') ':00']);
      delete(findobj(0,'Tag','uigetdate'))                       
   catch
      out = [];
      closereq
   end% try
   
   return
end% if

switch varargin{2}
   case 'months'
      h = findobj(gcbf,'Tag','months');
      months = get(h,'UserData');
      set(h,'String',months{mod(get(gcbo,'Value')-1,12)+1})
      set(findobj(gcbf,'Tag','ok'),'Enable','off')      
      %
   case 'year'
      set(findobj(gcbf,'Tag','year'),'String',get(gcbo,'Value'))
      set(findobj(gcbf,'Tag','ok'),'Enable','off')
      %
   case 'day'
      h= findobj(gcf,'Tag','day');
      set(h,'CData',[])

      set(varargin{1},'CData',geticon)
      set(findobj(gcbf,'Tag','cday'),'Userdata',get(varargin{1},'String'))
      set(findobj(gcbf,'Tag','ok'),'Enable','on')
      try ; uicontrol(h(3)) ; end% try
      return
      %
   case 'time'
      try
         if toc<0.1
            step = get(gcbo,'UserData');
            set(gcbo,'UserData',step+1)
            step = floor(step*sign(get(gcbo,'value'))/2);
         else
            set(gcbo,'UserData',1)
            step = sign(get(gcbo,'value'));
            set(gcbo,'value',0)
         end% if
         %
         handles.time = findobj(gcbf,'Tag','time');
         time = sum(sscanf(get(handles.time,'String'),'%d:%d').*[60;1]);
         time = time+step;
         if time<0
            time = 1439;
         elseif time>1439
            time = 0;
         end% if
         time = sprintf('%02.f:%02.f',floor(time/60),(time/60-floor(time/60))*60);
         set(handles.time,'String',time)
         %
         tic
         return
      catch
         tic
      end% try
      drawnow
      %
   case 'ok'
      uiresume
      return
      %
end% switch
setday(['1-' get(findobj(gcbf,'Tag','months'),'String') '-' ...
             get(findobj(gcbf,'Tag','year'),'String')])
  %
  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setday(datvec)
%-------------------------------------------------------------------------------
datvec = datevec(datvec);
datvec(3) = 1;
%
day = [7 1 2 3 4 5 6];
day = day(weekday(datestr(datvec)));
%
monthend = eomday(datvec(1),datvec(2));
%
ind = [zeros(1,42-monthend-day+1) monthend:-1:1 zeros(1,day-1)];
%
enable = repmat({'on'},42,1);
enable(ind==0) = {'off'};
%
count = strrep(strrep(cellstr(num2str(ind')),' 0',''),' ','');
%
h = findobj(gcf,'Tag','day');
set(h,{'String'},count,{'Enable'},enable,'backgroundcolor',[0.7 0.7 0.7],'CData',[])
set(h(ind~=0),'backgroundcolor',[.925 .922 .9002]);
set(h(ind~=0&repmat([1 1 0 0 0 0 0],1,6)),'backgroundcolor',[1 .8 .8])
  %
  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function icon = geticon
%-------------------------------------------------------------------------------
tmp = [0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ;
       0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 1 ; ...
       0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 1 1 ; ...
       1 1 0 0 0 0 1 1 1 1 1 0 0 0 1 1 1 ; ...
       1 1 0 0 0 0 0 1 1 0 0 0 0 1 1 1 1 ; ...
       1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 ; ...
       1 1 1 1 0 0 0 0 0 0 0 1 1 1 1 1 1 ; ...
       1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 ; ...
       1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 ; ...
       1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 ; ...
       1 1 0 0 0 0 0 1 1 0 0 0 0 0 1 1 1 ; ...
       0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 1 1 ; ...
       0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 1 1 ; ...
       0 0 0 1 1 1 1 1 1 1 1 1 0 0 0 0 1 ; ...
       1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 ];
tmp(tmp==1)=NaN;
tmp(tmp==0)=1;
icon(:,:,1) = tmp;
tmp(tmp==1)=0.25;
icon(:,:,2) = tmp;
tmp(tmp==.25)=0;
icon(:,:,3) = tmp;


