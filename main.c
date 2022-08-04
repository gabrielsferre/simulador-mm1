#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

enum TipoEvento {
    EVENTO_INVALIDO,
    EVENTO_CHEGADA_NA_FILA,
    EVENTO_SERVICO_COMPLETO,
};

enum Distribuicao {
    DISTRIBUICAO_INVALIDA,
    DISTRIBUICAO_CONSTANTE,
    DISTRIBUICAO_EXPONENCIAL,
    DISTRIBUICAO_DETERMINISTICA,
};

enum Disciplina {
    FILA_INVALIDA,
    FILA_FCFS,
    FILA_LCFS,
};

typedef struct Metricas_ {
    double somatorio_tempo_espera;
    double somatorio_tempo_espera_quadrado;
    double somatorio_media_tempo_espera;
    double somatorio_media_tempo_espera_quadrado;
    double somatorio_variancia_tempo_espera;
    double somatorio_variancia_tempo_espera_quadrado;
    double area_Nq_versus_tempo;
    double area_Nq_quadrado_versus_tempo;
    double somatorio_media_Nq;
    double somatorio_media_Nq_quadrado;
    double somatorio_variancia_Nq;
    double somatorio_variancia_Nq_quadrado;
    double ultima_vez_que_atualizou_area;
    double inicio_rodada;
    int coletas_clientes;
    int coletas_por_rodada;
    int rodada; // rodada zero é a fase transiente
    int numero_de_rodadas;
} Metricas;

typedef struct Evento_ {
    double tempo; // tempo em que o evento ocorreu
    enum TipoEvento tipo;
} Evento;

typedef struct ListaEventos_ {
    Evento array[1024];
    int tamanho_array;
    int tamanho_lista;
    int chegadas_sem_servico;
} ListaEventos;

typedef struct Cliente_ { // cliente da fila
    double chegada; // momento em que o cliente chegou na fila
    double comeco_servico;
} Cliente;

typedef struct Fila_ {
    Cliente array[1024];
    int tamanho_array;
    int inicio;
    int fim;
    int Nq;
    int servidor_ocupado;
    Cliente cliente_sendo_servido;
    Metricas metricas;
    enum Disciplina disciplina;
} Fila;

typedef struct Intervalo_ {
    double inf; // limite inferior do intervalo
    double sup; // limite superior do intervalo
} Intervalo;

typedef struct ParametrosDistribuicao_ {
    enum Distribuicao servidor;
    enum Distribuicao chegada;
    double parametro_servidor;
    double parametro_chegada;
} ParametrosDistribuicao;

typedef struct ParametrosSimulacao_ {
    enum Disciplina disciplina_fila;
    ParametrosDistribuicao distribuicao;
    int numero_de_rodadas;
    int coletas_por_rodada;
    int coletas_fase_transiente;
} ParametrosSimulacao;

double
calc_media_Nq(double rho)
{
    return (rho*rho)/(1-rho);
}

double
calc_variancia_Nq(double rho)
{
    return rho/((1-rho)*(1-rho)) - rho - rho*rho;
}

double
calc_media_tempo_espera(double rho)
{
    return rho/(1-rho);
}

double
calc_variancia_tempo_espera_FCFS(double rho)
{
    return (2*rho - rho*rho)/((1-rho)*(1-rho));
}

double
calc_variancia_tempo_espera_LCFS(double rho)
{
    double media = calc_media_tempo_espera(rho);
    return (calc_variancia_tempo_espera_FCFS(rho) + rho*media*media)/(1-rho);
}

double
precisao_intervalo(Intervalo i)
{
    return (i.sup - i.inf) / (i.sup + i.inf);
}

Intervalo
intervalo_chiquadrado(double variancia, int k)
{
    // alfa = 0.05 e 3199 graus de liberdade
    double quantil_1 = 3357.658; // probabilidade 1 - alfa/2
    double quantil_2 = 3044.13;  // probabilidade alfa/2
    double x = k*3199*variancia;
    Intervalo i = {.inf = x/quantil_1, .sup = x/quantil_2};
    return i;
}

Intervalo
intervalo_tstudent(double media, double variancia)
{
    // alfa=0.05 e 3199 graus de liberdade
    double quantil = 1.960706;
    double x = quantil*sqrt(variancia/3200.);
    Intervalo i = {.inf = media-x, .sup = media+x};
    return i;
}

double
estima_variancia(double soma_quadrados, double soma, int n)
{
    assert(n > 1);
    return (1./(n-1))*(soma_quadrados - soma*soma/n);
}

void
atualiza_metricas_fim_da_rodada(Metricas *m)
{
    // finaliza rodada
    if(m->coletas_clientes >= m->coletas_por_rodada) {
        assert(m->coletas_clientes == m->coletas_por_rodada);
        // atualiza metricas
        {
            double media_tempo_espera =
                m->somatorio_tempo_espera / m->coletas_clientes;
            m->somatorio_media_tempo_espera += media_tempo_espera;
            m->somatorio_media_tempo_espera_quadrado +=
                media_tempo_espera*media_tempo_espera;

            double variancia_tempo_espera =
                estima_variancia(m->somatorio_tempo_espera_quadrado,
                        m->somatorio_tempo_espera, m->coletas_clientes);
            m->somatorio_variancia_tempo_espera += variancia_tempo_espera;
            m->somatorio_variancia_tempo_espera_quadrado +=
                variancia_tempo_espera * variancia_tempo_espera;
            //printf("%i variancia: %f\n", m->rodada, variancia_tempo_espera);

            double dt = m->ultima_vez_que_atualizou_area - m->inicio_rodada; // TODO: checar se faz sentido fazer '(ultima_vez_que_atualizou_area - inicio_rodada)'
            double media_Nq = m->area_Nq_versus_tempo / dt;
            double media_Nq_quadrado = media_Nq*media_Nq;
            m->somatorio_media_Nq += media_Nq;
            m->somatorio_media_Nq_quadrado += media_Nq_quadrado;

            double variancia_Nq =
                m->area_Nq_quadrado_versus_tempo/dt - media_Nq_quadrado;
            m->somatorio_variancia_Nq += variancia_Nq;
            m->somatorio_variancia_Nq_quadrado += variancia_Nq*variancia_Nq;
        }
        // reinicia metricas da rodada
        {
            m->somatorio_tempo_espera = 0;
            m->somatorio_tempo_espera_quadrado = 0;
            m->area_Nq_versus_tempo = 0;
            m->area_Nq_quadrado_versus_tempo = 0;
            m->coletas_clientes = 0;
            m->inicio_rodada = m->ultima_vez_que_atualizou_area; // TODO: checar se isso faz sentido
        }
        m->rodada += 1;
    }
}

void
atualiza_metricas_tempo_espera(Cliente c, Metricas *m)
{
    double tempo_espera = c.comeco_servico - c.chegada;
    assert(tempo_espera >= 0);
    m->somatorio_tempo_espera += tempo_espera;
    m->somatorio_tempo_espera_quadrado +=
        tempo_espera*tempo_espera;
    m->coletas_clientes += 1;
}

void
atualiza_metricas_clientes_na_espera(double tempo, int Nq, Metricas *m)
{
    double dt = tempo - m->ultima_vez_que_atualizou_area;
    m->area_Nq_versus_tempo += Nq * dt;
    m->area_Nq_quadrado_versus_tempo += Nq*Nq*dt;
    m->ultima_vez_que_atualizou_area = tempo;
}

int
fila_proxima_posicao(int i, Fila const *fila)
{
    return (i+1 >= fila->tamanho_array) ? 0 : i+1; 
}

int
fila_posicao_anterior(int i, Fila const *fila)
{
    return (i-1 < 0) ? fila->tamanho_array - 1 : i-1; 
}

void
enfileira(Cliente c, Fila *fila)
{
    int proxima = fila_proxima_posicao(fila->fim, fila);
    assert(proxima != fila->inicio && "fila cheia");
    fila->array[fila->fim] = c;
    fila->fim = proxima;
}

Cliente
desenfileira(Fila *fila)
{
    assert(fila->inicio != fila->fim && "fila vazia");
    Cliente c = {0};
    if(fila->disciplina == FILA_FCFS) {
        c = fila->array[fila->inicio];
        fila->inicio = fila_proxima_posicao(fila->inicio, fila);
    }
    else if(fila->disciplina == FILA_LCFS) {
        int p = fila_posicao_anterior(fila->fim, fila);
        c = fila->array[p];
        fila->fim = p;
    }
    else {
        assert(0 && "disciplina invalida");
    }
    return c;
}

void
avanca_simulacao(Fila *fila, ListaEventos *eventos, Metricas *metricas)
{
    for(int i = 0;
            i < eventos->tamanho_lista &&
            metricas->rodada < metricas->numero_de_rodadas;
            ++i) {
        Evento e = eventos->array[i];
        atualiza_metricas_clientes_na_espera(e.tempo, fila->Nq, metricas);
        if(e.tipo == EVENTO_CHEGADA_NA_FILA) {
            Cliente c = {.chegada = e.tempo};
            if(fila->servidor_ocupado) {
                enfileira(c, fila);
                fila->Nq += 1;
            }
            else {
                c.comeco_servico = e.tempo;
                fila->cliente_sendo_servido = c;
                fila->servidor_ocupado = 1;
            }
        }
        else if(e.tipo == EVENTO_SERVICO_COMPLETO) {
            assert(fila->servidor_ocupado);
            Cliente cs = fila->cliente_sendo_servido;
            atualiza_metricas_tempo_espera(cs, metricas);
            atualiza_metricas_fim_da_rodada(metricas);
            bzero(&fila->cliente_sendo_servido, sizeof(Cliente));
            if(fila->Nq > 0) {
                Cliente c = desenfileira(fila);
                c.comeco_servico = e.tempo;
                fila->cliente_sendo_servido = c;
                fila->Nq -= 1;
            }
            else {
                fila->servidor_ocupado = 0;
            }

#if 0
            int iteracao = metricas->rodada*metricas->coletas_por_rodada +
                metricas->coletas_clientes;
            printf("%i Nq: %i\n", iteracao, fila->Nq);
#endif
        }
    }
}

void
insere(Evento e, ListaEventos *lista)
{
    assert(lista->tamanho_lista + 1 <= lista->tamanho_array);

    int i = lista->tamanho_lista;
    for(; i > 0; --i) {
        if(lista->array[i-1].tempo > e.tempo) {
            lista->array[i] = lista->array[i-1];
        }
        else {
            break;
        }
    }
    lista->array[i] = e;
    lista->tamanho_lista += 1;
}

double
random_01() {
    return random() / (double)((1ll << 31) - 1); // algo dentro do intervalo [0,1]
}

double
gera_amostra(enum Distribuicao distribuicao, double parametro)
{
    if(distribuicao == DISTRIBUICAO_CONSTANTE) {
        return parametro;
    }

    if(distribuicao == DISTRIBUICAO_DETERMINISTICA) {
        static int i = 0;
        i += 1;
        return (i%2 == 0) ? 4*parametro : parametro;
    }

    double u = random_01();
    
    if(distribuicao == DISTRIBUICAO_EXPONENCIAL) {
        double r = log(u);
        while(r == -HUGE_VAL) {
            u = random_01();
            r = log(u);
        }
        return -r/parametro;
    }
    else {
        assert(0 && "distribuicao invalida");
        return 0;
    }
}

void
print_info_tstudent(double valor_teorico,
        double somatorio, double somatorio_quadrados, int rodadas)
{
    double media = somatorio / rodadas;
    double variancia = estima_variancia(somatorio_quadrados, somatorio, rodadas);
    Intervalo intervalo = intervalo_tstudent(media, variancia);
    double precisao = precisao_intervalo(intervalo);
    int teorico_dentro_do_intervalo =
        valor_teorico >= intervalo.inf && valor_teorico <= intervalo.sup;

    char const *naosim[] = {"não", "sim"};
    char const *sp = "|   |   |-- ";
    printf("%sintervalo de confiança: [%f,%f]\n", sp, intervalo.inf, intervalo.sup);
    printf("%sprecisao do intervalo: %f\n", sp, precisao);
    printf("%scentro do intervalo: %f\n", sp, (intervalo.inf+intervalo.sup)/2);
    printf("%svalor teórico (analítico): %f\n", sp, valor_teorico);
    printf("%svalor teórico dentro do intervalo: %s\n", sp,
            naosim[teorico_dentro_do_intervalo]);
    printf("%svariancia do valor médio: %f\n", sp, variancia);
}

void
print_info_chiquadrado(double somatorio, double somatorio_quadrados,
        int rodadas, int coletas_por_rodada)
{
    double variancia = estima_variancia(somatorio_quadrados, somatorio, rodadas);
    Intervalo intervalo= intervalo_chiquadrado(variancia, coletas_por_rodada);
    double precisao = precisao_intervalo(intervalo);

    char const *sp = "|   |   |-- ";
    printf("%sintervalo de confiança: [%f,%f]\n", sp, intervalo.inf, intervalo.sup);
    printf("%sprecisao do intervalo: %f\n", sp, precisao);
    printf("%scentro do intervalo: %f\n", sp, (intervalo.inf+intervalo.sup)/2);
    printf("%svariancia do valor médio: %f\n", sp, variancia);
}

void
print_info(Metricas metricas, double rho, enum Disciplina disciplina_fila)
{
    // TODO: o cálculo da variancia usando chi quadrado está zaralhado, tenho q consertar isso
    assert(metricas.rodada == metricas.numero_de_rodadas);
    int rodadas = metricas.numero_de_rodadas;
    int coletas = metricas.coletas_por_rodada;

    char const *sp1 = "|-- ";
    printf("%sutilização usada: %.2f\n", sp1, rho);

    char const *sp2 = "|   |-- ";
    printf("%sMédia do tempo de espera:\n", sp2);
    print_info_tstudent(calc_media_tempo_espera(rho),
            metricas.somatorio_media_tempo_espera,
            metricas.somatorio_media_tempo_espera_quadrado, rodadas);
    printf("|   |\n");
    double variancia_teorica_tempo_espera = (disciplina_fila == FILA_FCFS) ?
        calc_variancia_tempo_espera_FCFS(rho) : calc_variancia_tempo_espera_LCFS(rho);
    printf("%sVariancia do tempo de espera usando t-student:\n", sp2);
    print_info_tstudent(variancia_teorica_tempo_espera,
            metricas.somatorio_variancia_tempo_espera,
            metricas.somatorio_variancia_tempo_espera_quadrado, rodadas);
    printf("|   |\n");
    printf("%sVariancia do tempo de espera usando chi quadrado:\n", sp2);
    print_info_chiquadrado(metricas.somatorio_media_tempo_espera,
            metricas.somatorio_media_tempo_espera_quadrado, rodadas, coletas);
    printf("|   |\n");

    printf("%sMédia de Nq:\n", sp2);
    print_info_tstudent(calc_media_Nq(rho),
            metricas.somatorio_media_Nq,
            metricas.somatorio_media_Nq_quadrado, rodadas);
    printf("|   |\n");
    printf("%sVariancia de Nq t-student:\n", sp2);
    print_info_tstudent(calc_variancia_Nq(rho),
            metricas.somatorio_variancia_Nq,
            metricas.somatorio_variancia_Nq_quadrado, rodadas);
    printf("|   |\n");
    printf("%sVariancia de Nq usando chi quadrado:\n", sp2);
    print_info_chiquadrado(metricas.somatorio_media_Nq,
            metricas.somatorio_media_Nq_quadrado, rodadas, coletas);
    printf("|   |\n");

    double media_espera = metricas.somatorio_media_tempo_espera / rodadas;
    double media_Nq = metricas.somatorio_media_Nq / rodadas;

    printf("%sestimativa da taxa de chegada usando Little: %.9f\n",
            sp2, media_Nq/media_espera);
    printf("%snumero de rodadas: %i\n", sp2, metricas.rodada);
    printf("%scoletas por rodada: %i\n", sp2, metricas.coletas_por_rodada);
    printf("%stempo considerado: %f\n", sp2, metricas.ultima_vez_que_atualizou_area);
}

double
loop_simulacao(ListaEventos *eventos, Fila *fila, Metricas *metricas,
        ParametrosDistribuicao dist_param, double tempo_inicial)
{
    assert(eventos->tamanho_lista == 1);
    for(double t_servico = tempo_inicial, t_chegada = tempo_inicial;
            metricas->rodada < metricas->numero_de_rodadas;) {

        assert(eventos->chegadas_sem_servico >= 1);
        if(eventos->chegadas_sem_servico == 1) {
            t_servico = t_chegada;
        }
        Evento servico = {
            .tempo = t_servico +
                gera_amostra(dist_param.servidor, dist_param.parametro_servidor),
            .tipo = EVENTO_SERVICO_COMPLETO,
        };
        insere(servico, eventos);
        eventos->chegadas_sem_servico -= 1;
        t_servico = servico.tempo;

        while(t_chegada <= t_servico) {
            Evento chegada = {0};
            chegada.tempo = t_chegada +
                gera_amostra(dist_param.chegada, dist_param.parametro_chegada);
            chegada.tipo = EVENTO_CHEGADA_NA_FILA;
            t_chegada = chegada.tempo;
            insere(chegada, eventos);
            eventos->chegadas_sem_servico += 1;
        }

        Evento ultima_chegada = eventos->array[eventos->tamanho_lista-1];
        assert(ultima_chegada.tipo == EVENTO_CHEGADA_NA_FILA);
        bzero(&eventos->array[eventos->tamanho_lista-1], sizeof(Evento));
        eventos->tamanho_lista -= 1;
        avanca_simulacao(fila, eventos, metricas);
        bzero(eventos->array, eventos->tamanho_lista*sizeof(Evento));
        eventos->tamanho_lista = 0;
        insere(ultima_chegada, eventos);
    }
    assert(eventos->tamanho_lista == 1);
    return eventos->array[0].tempo;
}

Metricas
gera_metricas(ParametrosSimulacao param)
{
    ListaEventos eventos = {0};
    Fila fila = {0};

    eventos.tamanho_array = sizeof(eventos.array)/sizeof(*eventos.array);
    fila.tamanho_array = sizeof(fila.array)/sizeof(*fila.array);
    fila.disciplina = param.disciplina_fila;

    Metricas metricas = {0};
    metricas.numero_de_rodadas = 1;
    metricas.coletas_por_rodada = param.coletas_fase_transiente;

    // gera primeira chegada
    {
        Evento chegada = {
            .tempo = 0,
            .tipo = EVENTO_CHEGADA_NA_FILA,
        };
        insere(chegada, &eventos);
        eventos.chegadas_sem_servico = 1;
    }
    // roda fase transiente
    double t = loop_simulacao(&eventos, &fila, &metricas, param.distribuicao, 0);

    bzero(&metricas, sizeof(metricas));
    metricas.numero_de_rodadas = param.numero_de_rodadas;
    metricas.coletas_por_rodada = param.coletas_por_rodada;
    loop_simulacao(&eventos, &fila, &metricas, param.distribuicao, t);

    return metricas;
}

int
main()
{
    ParametrosSimulacao param = {0};
    param.numero_de_rodadas = 3200;
    param.coletas_por_rodada = 50000;
    param.coletas_fase_transiente = 2000;
    param.distribuicao.servidor = DISTRIBUICAO_EXPONENCIAL;
    param.distribuicao.chegada = DISTRIBUICAO_EXPONENCIAL;
    param.distribuicao.parametro_servidor = 1;
    
    char const *disciplinas[] = {"FCFS", "LCFS"};
    for(int i = 0; i <= 1; ++i) {

        printf("Disciplina da fila: %s\n", disciplinas[i]);

        param.disciplina_fila = i == 0 ? FILA_FCFS : FILA_LCFS;

        double const lista_rho[] = {.2, .4, .6, .8, .9};
        int tamanho_lista_rho = sizeof(lista_rho)/sizeof(*lista_rho);
        for(int j = 0; j < tamanho_lista_rho; ++j) {
            double rho = lista_rho[j];
            param.distribuicao.parametro_chegada = rho;
            Metricas metricas = gera_metricas(param);
            printf("|\n");
            print_info(metricas, rho, param.disciplina_fila);
        }
        printf("\n");
    }
}
