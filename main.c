#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

enum TipoEvento {
    EVENTO_INVALIDO,
    EVENTO_CHEGADA_NA_FILA,
    EVENTO_SERVICO_COMPLETO,
};

enum Distribuicao {
    DISTRIBUICAO_INVALIDA,
    DISTRIBUICAO_CONSTANTE,
    DISTRIBUICAO_EXPONENCIAL,
};

enum Disciplina {
    FILA_INVALIDA,
    FILA_FCFS,
    FILA_LCFS,
};

typedef struct Metricas_ {
    double somatorio_tempo_espera_rodada;
    double somatorio_media_tempo_espera;
    double somatorio_media_tempo_espera_quadrado;
    double area_Nq_versus_tempo;
    double somatorio_media_Nq;
    double ultima_vez_que_atualizou_area;
    double inicio_rodada;
    int coletas_clientes;
    int coletas_por_rodada;
    int rodada;
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

void
atualiza_metricas_fim_da_rodada(Metricas *m)
{
    assert(m->coletas_clientes <= m->coletas_por_rodada);
    // finaliza rodada
    if(m->coletas_clientes >= m->coletas_por_rodada) {
        // atualiza metricas
        {
            double media_tempo_espera =
                m->somatorio_tempo_espera_rodada / m->coletas_clientes;
            m->somatorio_media_tempo_espera += media_tempo_espera;
            m->somatorio_media_tempo_espera_quadrado =
                media_tempo_espera*media_tempo_espera;

            double media_Nq = m->area_Nq_versus_tempo /
                (m->ultima_vez_que_atualizou_area - m->inicio_rodada); // TODO: checar se faz sentido fazer '(ultima_vez_que_atualizou_area - inicio_rodada)'
            m->somatorio_media_Nq += media_Nq;
        }
        // reinicia metricas da rodada
        {
            m->somatorio_tempo_espera_rodada = 0;
            m->area_Nq_versus_tempo = 0;
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
    m->somatorio_tempo_espera_rodada += tempo_espera;
    m->coletas_clientes += 1;
}

void
atualiza_metricas_clientes_na_espera(double tempo, int Nq, Metricas *m)
{
    double dt = tempo - m->ultima_vez_que_atualizou_area;
    m->area_Nq_versus_tempo += Nq * dt;
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
    }
}

int
main()
{
    ListaEventos *eventos = (ListaEventos*)calloc(1,sizeof(*eventos)); // fazendo isso pq 'eventos' vai ter q ser ponteiro depois
    Fila *fila = (Fila*)calloc(1,sizeof(*fila));
    Metricas *metricas = (Metricas*)calloc(1,sizeof(Metricas));

    eventos->tamanho_array = sizeof(eventos->array)/sizeof(*eventos->array);
    fila->tamanho_array = sizeof(fila->array)/sizeof(*fila->array);
    fila->disciplina = FILA_FCFS;
    metricas->coletas_por_rodada = 3500;
    metricas->numero_de_rodadas = 4000;

    // gera primeira chegada
    {
        Evento chegada = {
            .tempo = 0,
            .tipo = EVENTO_CHEGADA_NA_FILA,
        };
        insere(chegada, eventos);
        eventos->chegadas_sem_servico = 1;
    }

    for(double t_servico = 0, t_chegada = 0;
            metricas->rodada < metricas->numero_de_rodadas;) {

        assert(eventos->chegadas_sem_servico >= 1);
        if(eventos->chegadas_sem_servico == 1) {
            t_servico = t_chegada;
        }
        Evento servico = {
            .tempo = t_servico + gera_amostra(DISTRIBUICAO_EXPONENCIAL, 1),
            .tipo = EVENTO_SERVICO_COMPLETO,
        };
        insere(servico, eventos);
        eventos->chegadas_sem_servico -= 1;
        t_servico = servico.tempo;

        while(t_chegada <= t_servico) {
            Evento chegada = {0};
            chegada.tempo = t_chegada + gera_amostra(DISTRIBUICAO_EXPONENCIAL, .9),
            chegada.tipo = EVENTO_CHEGADA_NA_FILA,
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

    assert(metricas->rodada == metricas->numero_de_rodadas);
    double media_espera = metricas->somatorio_media_tempo_espera /
        metricas->numero_de_rodadas;
    double media_Nq = metricas->somatorio_media_Nq / 
        metricas->numero_de_rodadas;

    printf("media espera: %f\n", media_espera);
    printf("numero de rodadas: %i\n", metricas->rodada);
    printf("media Nq: %f\n", media_Nq);
    printf("tempo considerado: %f\n", metricas->ultima_vez_que_atualizou_area);
    printf("estimativa da taxa de chegada usando Little: %.9f\n", media_Nq/media_espera);
}
