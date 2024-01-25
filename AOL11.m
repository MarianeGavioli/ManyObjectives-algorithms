function[sol_final,Cost]=AOL1(Pop,params,ManyObj)

    tic;
       
    %% Behavioural parameters for this optimizer
    maxEvaluations = params.maxgen;
    NP = params.Np;        % Population size
    Dim = ManyObj.nVar;     % Variable number
    Bu = ManyObj.var_max;
    Bl = ManyObj.var_min;
    fun = ManyObj.fun;
    M = ManyObj.M;
    
    % Parameters
    % Separando a população entre: masc e fem
    % 70% serão leoas pride, 20% serão leoes pride e 10% leoes nomades
    Nleoas = 0.7*NP;
    Nleoes = 0.2*NP;
    Nleoes_nom = 0.1*NP;
    if Nleoes_nom <= 1
        Nleoes_nom = 2;
    end
    %%
    fitness = fun(Pop,M);
    [RANK] = FastNonDominatedSorting_Vectorized(fitness);
    [sortCrowd,sortRank,sortFit,sortPop] = crowdingDistances(RANK,fitness,Pop);
        
    % Inicializando população leão pride
    leoes = sortPop(1:Nleoes,:);
    Aptidao_leao = sortFit(1:Nleoes,:);
    % Inicializando população leoa pride
    leoas = sortPop(Nleoes+1:Nleoes+Nleoas,:);
    Aptidao_leoa = sortFit(Nleoes+1:Nleoes+Nleoas,:);
    % Inicializando população leão nomade
    leao_nom(1,:)= sortPop(Nleoas+Nleoes+1,:);
    Aptidao_leao_nom(1,:) = sortFit(Nleoas+Nleoes+1,:);
    %%    
    r12 = randi(Nleoes);
    fref = Aptidao_leao(r12,:);
    Lr = 0;
    Lrmax = 3;
    Sr = 0;
    Srmax = 3;
    %gcmax = round(0.3*Dim);
    gcmax = Dim;
    Agemax = 3;

    % ITERAÇÕES
       for it = 1:maxEvaluations
%--------------------------------------------------------------------------
            % Fertilização
%--------------------------------------------------------------------------
            % trazero lider como no mopso
            r12 = randi(Nleoes);
            d0 = dominates(fref,Aptidao_leao(r12,:));
            if d0 == 1
                Lr = Lr + 1;
            else
                Lr = 0;
                fref =  Aptidao_leao(r12,:);
            end
            if Sr <= Srmax
                uc = 0;
%                 gc = 0;
                for i1=1:gcmax
                    r11 = randi(Nleoas);
                    leoa_mais = leoas(r11,:);
                    r1 = rand();
                    r2 = rand();
                    k = randi(Dim);
                    delta = leoas(r11,k) + (0.1*r2 - 0.05)*(leoes(r12,k) - r1*leoas(r11,k));
                    leoa_mais(:,k) = min(Bu(1,k),max(Bl(1,k),delta));
                    Aptidao_leoa_mais = fun(leoa_mais,M);
                    d1 = dominates(Aptidao_leoa_mais,Aptidao_leoa(r11,:));
                    if d1 == 1
                        uc = uc + 1;
                        leoas(r11,:) = leoa_mais(:,:);
                        Aptidao_leoa(r11,:) = Aptidao_leoa_mais(:,:);
                    end
                end
                if uc == 0
                    Sr = Sr + 1;
                else
                    Sr = 0;
                end
            end
%--------------------------------------------------------------------------            
            % Acasalamento
%--------------------------------------------------------------------------
            % 1) Cruzamento
            j = 1;
            for i2=1:Nleoes
                c1 = randi(Dim);
                c2 = randi(Dim);
                ind = randi(Nleoas);
                cruzamento(j,:)   = [leoas(ind,1:c1) leoes(i2,c1+1:end)];
                cruzamento(j+1,:) = [leoas(ind,c1+1:end) leoes(i2,1:c1)];
                cruzamento(j+2,:) = [leoas(ind,1:c2) leoes(i2,c2+1:end)];
                cruzamento(j+3,:) = [leoas(ind,c2+1:end) leoes(i2,1:c2)];
                j = j + 4;
            end
                        
            % 2) Mutação
            mutacao = cruzamento;
            for k1 = 1:size(cruzamento,1)
                ind3 =  randi(Dim);
                mutacao(k1,ind3) = unifrnd(Bl(ind3),Bu(ind3));
            end
            
            % Avaliação dos filhotes:
            Cubs = [cruzamento; mutacao];
            clear cruzamento mutacao
            f_Cubs = fun(Cubs,M);
            [RANK_c] = FastNonDominatedSorting_Vectorized(f_Cubs);
            indxc = find(RANK_c == 1);
            cubs = Cubs(indxc,:);
            f_cubs = f_Cubs(indxc,:);
            Age = 0;
            defesa = 0;
            
%--------------------------------------------------------------------------            
            % Função de crescimento dos filhotes
%--------------------------------------------------------------------------
            while (Age < Agemax && defesa == 0)
                newMutation = normrnd(cubs,0.1);
                Flag2Bl = newMutation < Bl;
                [IMBl,~,~] = find(Flag2Bl == 1);
                newMutation(IMBl,:) = [];
                Flag2Bu = newMutation > Bu;
                [IMBu,~,~] = find(Flag2Bu == 1);
                newMutation(IMBu,:) = [];
                f_mutacao = fun(newMutation,M);
                pop_muta = [cubs; newMutation];
                fit_muta = [f_cubs; f_mutacao];
                clear Flag2Bu Flag2Bl
                [RANK_muta] = FastNonDominatedSorting_Vectorized(fit_muta);
                nsndm = find(RANK_muta == 1);
                cubs = pop_muta(nsndm,:);
                f_cubs = fit_muta(nsndm,:);
                RANK_fmuta = RANK_muta(nsndm);
                [sortCrowd_cubs,sortRank_fm,sortFit_cubs,sortcubs] = crowdingDistances1(RANK_fmuta,f_cubs,cubs);
                solrep = find(sortCrowd_cubs == 0);
                if numel(solrep) > 0
                    sortCrowd_cubs(solrep) = [];
                    sortRank_fm(solrep) = [];
                    sortFit_cubs(solrep,:) = [];
                    sortcubs(solrep,:) = [];
                end
                nFM = round(length(sortCrowd_cubs)/2);
                nFF = length(sortCrowd_cubs) - nFM;
%                 clear cubs f_cubs
%--------------------------------------------------------------------------               
        % Defesa de territorio
%--------------------------------------------------------------------------
                if Age == 0 && it == 1
                    num_nom = 1;
                else
                    num_nom = size(leao_nom,1);
                end
                r14 = randi(Nleoes,1,num_nom)'; %tem erro aqui
                d14 = dominates(Aptidao_leao_nom,Aptidao_leao(r14,:));
                ind14 = find(d14 == 1);
                if numel(ind14) > 0
                    Aptidao_leao(r14(ind14),:) = Aptidao_leao_nom(ind14,:);
                    leoes(r14(ind14),:) = leao_nom(ind14,:);
                    defesa = 1
                    Age = 0;    
                else
                    if Lr < Lrmax
                         for i = 2: Nleoes_nom
                             aux_leao_nom(i,:) = unifrnd(Bl,Bu);
                         end
                    else
                        for i = 2:Nleoes_nom
                            r13 = randi(Nleoes);
                            aux_leao_nom(i,:) = normrnd(leoes(r13,:),0.01);
                        end
                        Flag2Bl = aux_leao_nom < Bl;
                        [INBl,~,~] = find(Flag2Bl == 1);
                        aux_leao_nom(INBl,:) = [];
                        if size(aux_leao_nom,1) == 0
                            aux_leao_nom(1,:) = unifrnd(Bl,Bu);
                        end
                        Flag2Bu = aux_leao_nom > Bu;
                        [INBu,~,~] = find(Flag2Bu == 1);
                        aux_leao_nom(INBu,:) = [];
                        clear Flag2Bu Flag2Bl
                    end
                    f_nom = fun(aux_leao_nom,M);
                    [RANK_nom] = FastNonDominatedSorting_Vectorized(f_nom);
                    indxnom = find(RANK_nom == 1);
                    leao_nom = aux_leao_nom(indxnom,:);
                    Aptidao_leao_nom = f_nom(indxnom,:);
                    defesa = 0;
                    Age = Age + 1;
                end
            end
            
        
        %% Aquisição territorial
            if Age >= Agemax
                for ii=1:nFM
                    ri1 = randi(Nleoes);
                    ri2 = randi(nFM);
                    d2 = dominates(sortFit_cubs(ri2,:),Aptidao_leao(ri1,:));
                    if  d2 == 1
                        leoes(ri1,:) = sortcubs(ri2,:);
                        Aptidao_leao(ri1,:) = sortFit_cubs(ri2,:);
                    end
                end
                for j1=1:nFF
                    rj1 = randi(Nleoes);
                    rj2 = randi(nFF);
                    d3 = dominates(sortFit_cubs(rj2,:),Aptidao_leoa(rj1,:));
                    if  d3 == 1
                        leoas(rj1,:) = sortcubs(rj2,:);
                        Aptidao_leoa(rj1,:) = sortFit_cubs(rj2,:);
                        Sr = 0;
                    end
                end
            end
            % Update Best Cost
%             Flag2Bu = leoes < Bl;
%             Flag2Bl = leoes > Bu;
%             leoes = (leoes.*(~(Flag2Bu+Flag2Bl)))+Bu.*Flag2Bu+Bl.*Flag2Bl;
%             Flag4Bu = leoas < Bl;
%             Flag4Bl = leoas > Bu;
%             leoas =(leoas.*(~(Flag4Bu+Flag4Bl)))+Bu.*Flag4Bu+Bl.*Flag4Bl;
            Pop2 = [leoes; leoas; sortcubs; leao_nom];
            Aptidao2 = [Aptidao_leao; Aptidao_leoa; sortFit_cubs; Aptidao_leao_nom];
            [RANK_2] = FastNonDominatedSorting_Vectorized(Aptidao2);
            [sortCrowd_Pop2,sortRank_Pop2,sortFit_Pop2,sortPop2] = crowdingDistances1(RANK_2,Aptidao2,Pop2);
            solrep = find(sortCrowd_Pop2 == 0);
            if numel(solrep) > 0
                sortCrowd_Pop2(solrep) = [];
                sortRank_Pop2(solrep) = [];
                sortFit_Pop2(solrep,:) = [];
                sortPop2(solrep,:) = [];
            end
                
            leoes = sortPop2(1:Nleoes,:);
            Aptidao_leao = sortFit_Pop2(1:Nleoes,:);
            leoas = sortPop2(Nleoes+1:Nleoes+Nleoas,:);
            Aptidao_leoa = sortFit_Pop2(Nleoes+1:Nleoes+Nleoas,:);
            

%             BestCost(it,:) = Apt
%             bestPar (it,:) = leoes(1,:);
it
        end %fim for it
    
    ttempo = toc;
    
    Aptidao2 = [Aptidao_leao; Aptidao_leoa];
    
    sol_final = leoes(1,:);
    Cost = Aptidao2;
end
    
   
    
