%% ================================================================
%  Exemplo de Código LDPC (n=12, k=4) com Decodificação Sum-Product
%  - G e H construídos manualmente (código LDPC pequeno)
%  - Mapeamento BPSK em canal AWGN
%  - Decodificação por Belief Propagation (Sum-Product)
%
%  Referência teórica: Lathi - Modern Digital and Analog Communication
%  Systems (seção de códigos LDPC, iteração entre nós de verificação
%  e nós de variável).
%
%  Autor: Aleksander, Arthur, Gabriel, Isabella (adaptado com ajuda do ChatGPT)
% ================================================================

clear; clc; close all;

%% Parâmetros do código LDPC
k = 4;              % número de bits de informação
r = 8;              % número de bits de paridade
n = k + r;          % comprimento do código

% -----------------------------------------------------------------
% Matriz P (k x r) escolhida de forma esparsa
% G = [I_k  P]   (geradora sistemática)
% H = [P^T I_r]  (verificadora de paridade esparsa -> LDPC)
% -----------------------------------------------------------------
P = [ ...
    1 0 0 1 0 1 0 0;  % linha 1
    0 1 0 1 1 0 0 0;  % linha 2
    0 0 1 0 1 0 1 0;  % linha 3
    1 0 0 0 0 1 0 1]; % linha 4

G = [eye(k) P];           % (4 x 12)
H = [P.' eye(r)];         % (8 x 12)

% Verificação rápida: H * G^T deve ser zero em GF(2)
teste = mod(H * G.', 2);
if any(teste(:) ~= 0)
    warning('H*G^T != 0 (algo errado na construção do código)');
else
    fprintf('OK: H*G^T = 0 (código consistente em GF(2)).\n');
end

%% Parâmetros da simulação
EbN0dB   = 3;       % Eb/N0 em dB (pode fazer sweep depois)
Nblocks  = 2000;    % número de palavras de código transmitidas
maxIter  = 20;      % iterações máximas do decodificador

% Relação entre Eb/N0 e variância do ruído
R        = k/n;                         % taxa do código
EbN0     = 10^(EbN0dB/10);
EsN0     = EbN0 * R;                    % Es/N0 (BPSK, Es=1)
sigma2   = 1/(2*EsN0);                  % variância do ruído
sigma    = sqrt(sigma2);

nBitErr  = 0;
nBitTot  = Nblocks * k;

%% Pré-processamento: listas de vizinhança (Tanner graph)
[M, N] = size(H);   % M linhas de checagem, N = n variáveis

% Para cada check node m, lista de variáveis conectadas
CN = cell(M,1);
for m = 1:M
    CN{m} = find(H(m,:) == 1);
end

% Para cada variable node n, lista de checks conectados
VN = cell(N,1);
for v = 1:N
    VN{v} = find(H(:,v) == 1).';
end

%% Loop de simulação
for blk = 1:Nblocks

    % --------------------------------------------------------------
    % 1) Geração de bits de informação (mensagem)
    % --------------------------------------------------------------
    u = randi([0 1], 1, k);         % mensagem binária

    % --------------------------------------------------------------
    % 2) Codificação LDPC: c = u * G (em GF(2))
    % --------------------------------------------------------------
    c = mod(u * G, 2);              % palavra de código de tamanho n

    % --------------------------------------------------------------
    % 3) Mapeamento BPSK: 0 -> +1, 1 -> -1
    % --------------------------------------------------------------
    x = 1 - 2*c;                    % vetor de símbolos BPSK (±1)

    % --------------------------------------------------------------
    % 4) Canal AWGN
    % --------------------------------------------------------------
    y = x + sigma * randn(size(x)); % recebido no canal

    % LLR inicial de cada bit (a priori)
    Lch = 2*y / sigma2;             % LLR do canal (log-likelihood ratio)

    % --------------------------------------------------------------
    % 5) Decodificador LDPC - Sum-Product (Belief Propagation)
    % --------------------------------------------------------------

    % Mensagens dos nós de variável para os nós de checagem: q(m,n)
    q = zeros(M,N);
    % Inicialmente, todos recebem somente a LLR do canal
    for m = 1:M
        for v = CN{m}
            q(m,v) = Lch(v);
        end
    end

    % Mensagens dos nós de checagem para os nós de variável: r(m,n)
    r_msg = zeros(M,N);

    % Iterações
    for it = 1:maxIter

        % ----- Atualização dos CHECK NODES (r_msg) -----
        for m = 1:M
            vlist = CN{m};
            for v = vlist
                % produto dos tanh das mensagens q/2 (exceto v)
                idx_others = vlist(vlist ~= v);
                prod_tanh = 1;
                for v2 = idx_others
                    prod_tanh = prod_tanh * tanh(q(m,v2)/2);
                end
                % mensagem do nó de verificação para o nó de variável
                r_msg(m,v) = 2 * atanh(prod_tanh + 1e-12);  % +eps p/ estabilidade numérica
            end
        end

        % ----- Atualização dos VARIABLE NODES (q) -----
        for v = 1:N
            mlist = VN{v};
            for m = mlist
                % soma das mensagens de todos os outros check nodes + canal
                idx_others = mlist(mlist ~= m);
                sum_r = sum(r_msg(idx_others, v));
                q(m,v) = Lch(v) + sum_r;
            end
        end

        % ----- Cálculo dos LLR "a posteriori" e decisão dura -----
        LLR_ap = zeros(1,N);
        for v = 1:N
            mlist = VN{v};
            LLR_ap(v) = Lch(v) + sum(r_msg(mlist, v));
        end

        c_hat = (LLR_ap < 0);   % decisão hard (0/1)

        % Verificação de síndrome: se H*c_hat^T = 0, sai mais cedo
        syndrome = mod(H * c_hat.', 2);
        if all(syndrome == 0)
            break;  % convergiu antes do número máximo de iterações
        end

    end % fim das iterações

    % --------------------------------------------------------------
    % 6) Contagem de erros somente nos bits de informação
    %    (primeiros k bits se G = [I_k P])
    % --------------------------------------------------------------
    u_hat = c_hat(1:k);
    nBitErr = nBitErr + sum(u ~= u_hat);

end % fim do loop de blocos

BER = nBitErr / nBitTot;

fprintf('\n===== RESULTADOS LDPC (n=%d, k=%d, R=%.3f) =====\n', n, k, R);
fprintf('Eb/N0 = %.1f dB,  Blocos = %d\n', EbN0dB, Nblocks);
fprintf('BER simulada = %.3e\n', BER);
fprintf('Iterações máx. = %d (com parada antecipada se síndrome = 0)\n', maxIter);
