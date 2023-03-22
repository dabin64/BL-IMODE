classdef smd5mp1
    properties
        p = 1;
        q = 2;
        r = 1;
        n_lvar;
        n_uvar;
        xu_bl;
        xu_bu;
        xl_bl;
        xl_bu;
        name;
        fu_prime;
        fl_prime ;
        xu_prime;
        xl_prime;
    end
    methods
        function obj = smd5mp1(p, q, r)
            if nargin == 3
                obj.p = p;
                obj.q = q;
                obj.r = r;
            end
            obj.name = 'SMD5mp1';
            
            % level variables
            obj.n_lvar = obj.q + obj.r;
            obj.n_uvar = obj.p + obj.r;
            
            % bounds
            %init bound upper level
            xu_bl_1 = ones(1, obj.p)* (-5.0);
            xu_bu_1 = ones(1, obj.p) * 10.0;
            xu_bl_2 = ones(1, obj.r)* (-5.0);
            xu_bu_2 = ones(1, obj.r) * 10.0;
            obj.xu_bl = [xu_bl_1, xu_bl_2];
            obj.xu_bu = [xu_bu_1, xu_bu_2];
            
            % init bound lower level
            xl_bl_1 = ones(1, obj.q) * (-5.0);
            xl_bu_1 = ones(1, obj.q) * 10.0;
          
            % xl_bl_2 = ones(1, obj.r) * (-5.0);
            % xl_bu_2 = ones(1, obj.r) * 10.0;
            
            xl_bl_2 = ones(1, obj.r) * 0;
            xl_bu_2 = ones(1, obj.r) * 10.0;
            
            
            obj.xl_bl = [xl_bl_1, xl_bl_2];
            obj.xl_bu = [xl_bu_1, xl_bu_2];
            
            obj.fu_prime = 0;
            xu_prime = zeros(1, obj.n_uvar);
            xl_prime = obj.get_xlprime(xu_prime);
            obj.fl_prime = obj.evaluate_l(xu_prime, xl_prime);
            
            
        end
        
        function [f, c] = evaluate_u(obj, xu, xl)
            xu1 = xu(:, 1 : obj.p);
            xu2 = xu(:, obj.p + 1: obj.p + obj.r);
            
            xl1 = xl(:, 1 : obj.q);
            xl2 = xl(:, obj.q+1 : obj.q+obj.r);
            
            term2 = 0;
            for i=1:obj.q-1
                term2 = term2 + (xl1(:, i+1) - xl1(:, i).^2).^2 + (xl1(:, i) - 1).^2;
            end
            
            xu_prime2 = zeros(1,  obj.r);            
            reg = 10;           % absolute value of  lower bound of xl2        
            t = 0.1 * tan(pi/2 - 0.01 - (abs(xu(:, obj.p+1 : end) )- xu_prime2) * pi/(2 * reg));
            
            %-obj
            f = sum((xu1).^2, 2) ...
                - term2 ...
                + sum((xu2).^2, 2) ...
                -  sum( (xl2 - t).^2, 2); 
                % - sum((abs(xu2) - xl2.^2).^2, 2);
 
            %-cie
            c = [];
        end
        
        function [f, c] = evaluate_l(obj, xu, xl)
            xu1 = xu(:, 1 : obj.p);
            xu2 = xu(:, obj.p + 1: obj.p + obj.r);
            
            xl1 = xl(:, 1 : obj.q);
            xl2 = xl(:, obj.q+1 : obj.q+obj.r);
            %-obj
            term2 = 0;
            for i=1:obj.q-1
                term2 = term2 + (xl1(:, i+1) - xl1(:, i).^2).^2 + (xl1(:, i) - 1).^2;
            end
            
            xu_prime2 = zeros(1,  obj.r);            
            reg = 10; % absolute value of  lower bound of xl2        
            t = 0.1 * tan(pi/2 - 0.01 - (abs(xu(:, obj.p+1 : end) )- xu_prime2) * pi/(2 * reg));
            
            f = sum((xu1).^2, 2) ...
                + term2 ...
                +  sum( (xl2 - t).^2, 2); 
                % + sum((abs(xu2) - xl2.^2).^2, 2);
            
             % grafting moving peak problem
%             fmp = [];
%             xl_center= obj.get_xlprime(xu);         
%             % xl_center = [0, 0.7];
%             
%             nx    = size(xu, 1);
%             w     = ones(5, 1) * 0.02;
%             h      = ones(5, 1) * 150;
%             h(1)  = 200;
% 
%             v = 0.3 .*  (obj.xl_bu(2) - obj.xl_bl(2));      % 1/10 of variable range
%             v = [0, v];
%             delta = [100/(obj.xl_bu(1) - obj.xl_bl(1)), 100/(obj.xl_bu(2) - obj.xl_bl(2))];
%             forward =  xl_center + v;
%             backward = xl_center - v;
%             
%              v2 = 0.3 .*  (obj.xl_bu(1) - obj.xl_bl(1));
%              v2 = [v2, 0];
%             left = xl_center + v2;
%             right =  xl_center - v2;
%             mp_centers = {xl_center, forward, backward, left, right};
%            
%             fmp = [];
%             for j = 1:length(mp_centers)
%                 fj = h(j) ./ (1 + w(j)  .* sum(((xl - mp_centers{j}) .*  delta).^2,  2)); 
%                 fmp = [fmp, fj];
%             end
%             
%             fmp = max(fmp, [], 2);
            fmp = mp_module(obj, xu , xl);
            f = f - fmp;
                       
            %-cie
            c = [];
            
            
        end
        
        function xl_prime = get_xlprime(obj, xu)
            n = size(xu, 1);
            xl1 = ones(n, obj.q);
            
            xu_prime2 = zeros(1, obj.r);            
            reg = 10; % absolute value of  lower bound of xl2        
            xl2 = 0.1 * tan(pi/2 - 0.01 - (abs(xu(:, obj.p+1 : end)) - xu_prime2) * pi/(2 * reg));
            
           
          
            xl = [xl1, xl2];
            xl_prime = xl;
        end
        
        function xl_local = get_otherlocal(obj, xu)
            xl_local = mp_othercenters(obj, xu);
        end
    end
end
