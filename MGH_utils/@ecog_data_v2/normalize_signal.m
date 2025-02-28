function  normalize_signal(obj,args)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments
    obj ecog_data_v2
    args.normtype = 'z-score'
end

if(~isfield(obj.stats,'normMetrics'))
    obj.extract_normalization_metrics();
end

normMetrics = obj.stats.normMetrics;

switch args.normtype
    case 'z-score'
        obj.elec_data = (obj.elec_data-normMetrics.normFactor(:,1))./normMetrics.normFactor(:,2);
    case 'mean-sub'
        obj.elec_data = (obj.elec_data-normMetrics.normFactor(:,1));
    case 'perc-change'
        obj.elec_data = (obj.elec_data-normMetrics.normFactor(:,1))./normMetrics.normFactor(:,1);
    case 'ratio'
        obj.elec_data = (obj.elec_data./normMetrics.normFactor(:,1));
    case 'log-ratio'
        obj.elec_data = 10.*log10(obj.elec_data./normMetrics.normFactor(:,1));
    case 'norm'
        obj.elec_data = (obj.elec_data-normMetrics.normFactor(:,1))./(obj.elec_data+normMetrics.normFactor(:,1));
end

if(~isempty(obj.bip_elec_data))
    switch args.normtype
        case 'z-score'
            obj.bip_elec_data = (obj.bip_elec_data-normMetrics.normFactor_bip(:,1))./normMetrics.normFactor_bip(:,2);
        case 'mean-sub'
            obj.bip_elec_data = (obj.bip_elec_data-normMetrics.normFactor_bip(:,1));
        case 'perc-change'
            obj.bip_elec_data = (obj.bip_elec_data-normMetrics.normFactor_bip(:,1))./normMetrics.normFactor_bip(:,1);
        case 'ratio'
            obj.bip_elec_data = (obj.bip_elec_data./normMetrics.normFactor_bip(:,1));
        case 'log-ratio'
            obj.bip_elec_data = 10.*log10(obj.bip_elec_data./normMetrics.normFactor_bip(:,1));
        case 'norm'
            obj.bip_elec_data = (obj.bip_elec_data-normMetrics.normFactor_bip(:,1))./(obj.bip_elec_data+normMetrics.normFactor_bip(:,1));
    end
end

end

