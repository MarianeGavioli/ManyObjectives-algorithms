function[sol_final,Cost]=AOL12(Pop,params,ManyObj)

    tic;
       
    % Behavioural parameters for this optimizer
    maxEvaluations = params.maxgen;
    NP = params.Np;        % Population size
    Dim = ManyObj.nVar;     % Variable number
    Bu = ManyObj.var_max;
    Bl = ManyObj.var_min;
    fun = ManyObj.fun;
    M = ManyObj.M;
    
    rep.ind = [];
    rep.Apt = [];
    rep_m.ind = [];
    rep_m.Apt = [];
    rep_f.ind = [];
    rep_f.Apt = [];
    
    % Parameters
    % Separando a população entre: masc e fem
    % 70% serão leoas pride, 20% serão leoes pride e 10% leoes nomades
    Nleoas = 0.7*NP;
    Nleoes = 0.2*NP;
    Nleoes_nom = 0.1*NP;
    if Nleoes_nom <= 1
        Nleoes_nom = 2;
    end
    
    fitness = fun(Pop,M);

    % Inicializando população leão pride
    leoes = Pop(1:Nleoes,:);
    Aptidao_leao = fitness(1:Nleoes,:);
    % Inicializando população leoas pride
    leoas = Pop(Nleoes+1:Nleoes+Nleoas,:);
    Aptidao_leoa = fitness(Nleoes+1:Nleoes+Nleoas,:);
    % Inicializando população leão nomade
    leao_nom(1,:)= Pop(NP,:);
    Aptidao_leao_nom(1,:) = fitness(NP,:);
    
    gcmax = Dim;
    Agemax = 3;

    % ITERAÇÕES
    %----------------------------------------------------------------------
    for it = 1:maxEvaluations
        % Fertilização
        %------------------------------------------------------------------
        for i1=1:gcmax
            r11 = randi(Nleoas);
            r12 = randi(Nleoes);
            leoa_mais = leoas(r11,:);
            r1 = rand();
            r2 = rand();
            k = randi(Dim);
            delta = leoas(r11,k) + (0.1*r2 - 0.05)*(leoes(r12,k) - r1*leoas(r11,k));
            leoa_mais(:,k) = min(Bu(1,k),max(Bl(1,k),delta));
            Aptidao_leoa_mais = fun(leoa_mais,M);
            d1 = dominates(Aptidao_leoa_mais,Aptidao_leoa(r11,:));
            d2 = dominates(Aptidao_leoa(r11,:),Aptidao_leoa_mais);
            if d1 == 1
                leoas(r11,:) = leoa_mais(:,:);
                Aptidao_leoa(r11,:) = Aptidao_leoa_mais(:,:);
            else
                if d2 == 0
                    rep_f.ind = [rep_f.ind; leoa_mais(:,:)];
                    rep_f.Apt = [rep_f.Apt; Aptidao_leoa_mais(:,:)];
                end
            end
        end
%--------------------------------------------------------------------------
        % Acasalamento
        %------------------------------------------------------------------
        % 1) Cruzamento
        j = 1;
        for i2=1:Nleoes
            c1 = randi([0 1],Dim,1);
            c2 = ~c1;
            
            c1 = randi(Dim);
            c2 = randi(Dim);
            ind(i2) = randi(Nleoas);
            cruzamento(j,:)   = [leoas(ind(i2),1:c1) leoes(i2,c1+1:end)];
            cruzamento(j+1,:) = [leoas(ind(i2),c1+1:end) leoes(i2,1:c1)];
            cruzamento(j+2,:) = [leoas(ind(i2),1:c2) leoes(i2,c2+1:end)];
            cruzamento(j+3,:) = [leoas(ind(i2),c2+1:end) leoes(i2,1:c2)];
            j = j + 4;
        end                
        % 2) Mutação
        mutacao = cruzamento;
        for k1 = 1:size(mutacao,1)
            ind3 =  randi(Dim);
            mutacao(k1,ind3) = unifrnd(Bl(ind3),Bu(ind3));
        end    
        % 3) Avaliação dos filhotes:
        Cubs = [cruzamento; mutacao];
        f_Cubs = fun(Cubs,M);
        nn = size(cruzamento,1);
        clear cruzamento mutacao
        j = 1;
        for ic=1:Nleoes
            cubs_m(j,:) =     Cubs(j,:);
            f_cubs_m(j,:) =   f_Cubs(j,:);
            cubs_m(j+1,:) =   Cubs(j+1,:);
            f_cubs_m(j+1,:) = f_Cubs(j+1,:);
            cubs_m(j+2,:) =   Cubs(nn+j,:);
            f_cubs_m(j+2,:) = f_Cubs(nn+j,:);
            cubs_m(j+3,:) =   Cubs(nn+j+1,:);
            f_cubs_m(j+3,:) = f_Cubs(nn+j+1,:);
            cubs_f(j,:) =     Cubs(j+2,:);
            f_cubs_f(j,:) =   f_Cubs(j+2,:);
            cubs_f(j+1,:) =   Cubs(j+3,:);
            f_cubs_f(j+1,:) = f_Cubs(j+3,:);
            cubs_f(j+2,:) =   Cubs(nn+j+2,:);
            f_cubs_f(j+2,:) = f_Cubs(nn+j+2,:);
            cubs_f(j+3,:) =   Cubs(nn+j+3,:);
            f_cubs_f(j+3,:) = f_Cubs(nn+j+3,:);
            j = j + 4;
        end           
        Age = 0;
        defesa = 0;
        cubs = [cubs_m; cubs_f];
        f_cubs = [f_cubs_m; f_cubs_f];
%--------------------------------------------------------------------------
        % Defesa de territorio
        %------------------------------------------------------------------
        % Função de crescimento dos filhotes
        while (Age < Agemax && defesa == 0)
            Age = Age + 1;
            newMutation = normrnd(cubs,0.1);
            newMutation = min(max(newMutation,Bl),Bu);
            f_mutacao = fun(newMutation,M);
            for inm=1:nn
                d11 = dominates(f_mutacao(inm,:),f_cubs(inm,:));
                d12 = dominates(f_cubs(inm,:),f_mutacao(inm,:));
                if d11 == 1
                    cubs(inm,:) = newMutation(inm,:);
                    f_cubs(inm,:) = f_mutacao(inm,:);
                else
                    if d12 == 0
                        rep_m.ind = [rep_m.ind; newMutation(inm,:)];
                        rep_m.Apt = [rep_m.Apt; f_mutacao(inm,:)];
                    end
                end
            end
            for inm=nn+1:(2*nn)
                d13 = dominates(f_mutacao(inm,:),f_cubs(inm,:));
                d14 = dominates(f_cubs(inm,:),f_mutacao(inm,:));
                if d13 == 1
                    cubs(inm,:) = newMutation(inm,:);
                    f_cubs(inm,:) = f_mutacao(inm,:);
                else
                    if d14 == 0
                        rep_f.ind = [rep_f.ind; newMutation(inm,:)];
                        rep_f.Apt = [rep_f.Apt; f_mutacao(inm,:)];
                    end
                end
            end
            % Ataque dos leoes nomades
            if it == 1
                num_nom = 1;
            else
                num_nom = size(leao_nom,1);
            end
            for inm = 1:num_nom
                r14 = randi(Nleoes);
                d15 = dominates(Aptidao_leao_nom(inm,:),Aptidao_leao(r14,:));
                d16 = dominates(Aptidao_leao(r14,:),Aptidao_leao_nom(inm,:));
                if d15 == 1
                    Aptidao_leao(r14,:) = Aptidao_leao_nom(inm,:);
                    leoes(r14,:) = leao_nom(inm,:);
                    defesa = 1;
                    Age = 0;
                else
                    if d16 == 0
                        rep_m.ind = [rep_m.ind; leao_nom(inm,:)];
                        rep_m.Apt = [rep_m.Apt; Aptidao_leao_nom(inm,:)];
                        defesa = 1;
                        Age = 0;
                    end
                end
            end
            if defesa == 0 && Age < Agemax
                for i = 1:Nleoes_nom
                    aux_leao_nom(i,:) = unifrnd(Bl,Bu);
                    r13 = randi(Nleoes);
                    aux_leao_nom1(i,:) = normrnd(leoes(r13,:),0.1);
                end
                aux_leao_nom=[aux_leao_nom; aux_leao_nom1];
                aux_leao_nom = min(max(aux_leao_nom,Bl),Bu);
                f_nom = fun(aux_leao_nom,M);
                [RANK_nom] = FastNonDominatedSorting_Vectorized(f_nom);
                indxnom = find(RANK_nom == 1);
                leao_nom = aux_leao_nom(indxnom,:);
                Aptidao_leao_nom = f_nom(indxnom,:);
                clear aux_leao_nom f_nom
                defesa = 0;
            end
        end
%--------------------------------------------------------------------------
        % Aquisição territorial
        %------------------------------------------------------------------
        if Age >= Agemax
            j = 1;
            for iat=1:Nleoes
                d21 = dominates(f_cubs_m(j,:),Aptidao_leao(iat,:));
                d22 = dominates(Aptidao_leao(iat,:),f_cubs_m(j,:));
                if d21 == 1
                    leoes(iat,:) = cubs_m(j,:);
                    Aptidao_leao(iat,:) = f_cubs_m(j,:);
                else
                    if d22 == 0
                        rep_m.ind = [rep_m.ind; cubs_m(j,:)];
                        rep_m.Apt = [rep_m.Apt; f_cubs_m(j,:)];
                    end
                end
                d31 = dominates(f_cubs_m(j+1,:),Aptidao_leao(iat,:));
                d32 = dominates(Aptidao_leao(iat,:),f_cubs_m(j+1,:));
                if d31 == 1
                    leoes(iat,:) = cubs_m(j+1,:);
                    Aptidao_leao(iat,:) = f_cubs_m(j+1,:);
                else
                    if d32 == 0
                        rep_m.ind = [rep_m.ind; cubs_m(j+1,:)];
                        rep_m.Apt = [rep_m.Apt; f_cubs_m(j+1,:)];
                    end
                end
                d41 = dominates(f_cubs_m(j+2,:),Aptidao_leao(iat,:));
                d42 = dominates(Aptidao_leao(iat,:),f_cubs_m(j+2,:));
                if d41 == 1
                    leoes(iat,:) = cubs_m(j+2,:);
                    Aptidao_leao(iat,:) = f_cubs_m(j+2,:);
                else
                    if d42 == 0
                        rep_m.ind = [rep_m.ind; cubs_m(j+2,:)];
                        rep_m.Apt = [rep_m.Apt; f_cubs_m(j+2,:)];
                    end
                end
                d51 = dominates(f_cubs_m(j+3,:),Aptidao_leao(iat,:));
                d52 = dominates(Aptidao_leao(iat,:),f_cubs_m(j+3,:));
                if d51 == 1
                    leoes(iat,:) = cubs_m(j+3,:);
                    Aptidao_leao(iat,:) = f_cubs_m(j+3,:);
                else
                    if d52 == 0
                        rep_m.ind = [rep_m.ind; cubs_m(j+3,:)];
                        rep_m.Apt = [rep_m.Apt; f_cubs_m(j+3,:)];
                    end
                end
                j = j+4;
            end
            k = 1;
            for iat=1:Nleoas
                d61 = dominates(f_cubs_f(k,:),Aptidao_leoa(ind(iat),:));
                d62 = dominates(Aptidao_leoa(iat,:),f_cubs_f(k,:));
                if d61 == 1
                    leoas(ind(iat),:) = cubs_f(k,:);
                    Aptidao_leoa(ind(iat),:) = f_cubs_f(k,:);
                else
                    if d62 == 0
                        rep_f.ind = [rep_f.ind; cubs_f(k,:)];
                        rep_f.Apt = [rep_f.Apt; f_cubs_f(k,:)];
                    end
                end
                d71 = dominates(f_cubs_f(k+1,:),Aptidao_leoa(ind(iat),:));
                d72 = dominates(Aptidao_leoa(ind(iat),:),f_cubs_f(k+1,:));
                if d71 == 1
                    leoas(ind(iat),:) = cubs_f(k+1,:);
                    Aptidao_leoa(ind(iat),:) = f_cubs_f(k+1,:);
                else
                    if d72 == 0
                        rep_f.ind = [rep_f.ind; cubs_f(k+1,:)];
                        rep_f.Apt = [rep_f.Apt; f_cubs_f(k+1,:)];
                    end
                end
                d81 = dominates(f_cubs_f(k+2,:),Aptidao_leoa(ind(iat),:));
                d82 = dominates(Aptidao_leoa(ind(iat),:),f_cubs_f(k+2,:));
                if d81 == 1
                    leoas(ind(iat),:) = cubs_f(k+2,:);
                    Aptidao_leoa(ind(iat),:) = f_cubs_f(k+2,:);
                else
                    if d82 == 0
                        rep_f.ind = [rep_f.ind; cubs_f(k+2,:)];
                        rep_f.Apt = [rep_f.Apt; f_cubs_f(k+2,:)];
                    end
                end
                d91 = dominates(f_cubs_f(k+3,:),Aptidao_leoa(ind(iat),:));
                d92 = dominates(Aptidao_leoa(ind(iat),:),f_cubs_f(k+3,:));
                if d91 == 1
                    leoas(ind(iat),:) = cubs_f(k+3,:);
                    Aptidao_leoa(ind(iat),:) = f_cubs_f(k+3,:);
                else
                    if d92 == 0
                        rep_f.ind = [rep_f.ind; cubs_f(k+3,:)];
                        rep_f.Apt = [rep_f.Apt; f_cubs_f(k+3,:)];
                    end
                end
                k = k+4;
            end
        end
%--------------------------------------------------------------------------            
        aux.ind = [rep.ind; leoes; leoas; rep_m.ind; rep_f.ind];
        aux.Apt = [rep.Apt; Aptidao_leao; Aptidao_leoa; rep_m.Apt; rep_f.Apt];
        [aux] = popnd(aux);
        rep_m.ind(:,:) = [];
        rep_f.ind(:,:) = [];
        rep_m.Apt(:,:) = [];
        rep_f.Apt(:,:) = [];
        [aux_Rank] = FastNonDominatedSorting_Vectorized(aux.Apt);
        [sortCrowd,sortRank,sortFit,sortPop] = crowdingDistances(aux_Rank,aux.Apt,aux.ind);
        indxf1 = find(sortRank == 1);
        rep.ind = sortPop(indxf1,:);
        rep.Apt = sortFit(indxf1,:);
        [rep] = popnd1(rep,NP);
        Nsol = size(rep.ind,1)
        leoes = sortPop(1:Nleoes,:);
        Aptidao_leao = sortFit(1:Nleoes,:);
        leoas = sortPop(Nleoes+1:Nleoes+Nleoas,:);
        Aptidao_leoa = sortFit(Nleoes+1:Nleoes+Nleoas,:);
it            
    end %fim for it
%--------------------------------------------------------------------------
    [rep] = popnd(rep);
    ttempo = toc;
    
    Aptidao2 = [Aptidao_leao; Aptidao_leoa];
    
    sol_final = leoes(1,:);
    Cost = Aptidao2;
end
    
   
    
