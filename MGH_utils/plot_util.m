classdef plot_util < handle
    %PLOT_UTIL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        
    %% method for constructing the object 
        function obj = plot_util(ecog_data)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.ecog_data = ecog_data;
        end    
        function obj = plot_pial(obj,inputArg2)
            %PLOT_UTIL Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

