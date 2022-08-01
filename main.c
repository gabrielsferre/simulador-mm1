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
    double somatorio_tempo_espera;
    double somatorio_tempo_espera_ao_quadrado;
    int coletas_clientes;
    double area_Nq_versus_tempo;
    double ultima_vez_que_atualizou_area;
} Metricas;

typedef struct Evento_ {
    double tempo; // tempo em que o evento ocorreu
    enum TipoEvento tipo;
} Evento;

typedef struct ListaEventos_ {
    Evento array[1024];
    int tamanho_array;
    int tamanho_lista;
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
atualiza_metricas_tempo_espera(Cliente c, Metricas *m)
{
    double tempo_espera = c.comeco_servico - c.chegada;
    assert(tempo_espera >= 0);
    m->somatorio_tempo_espera += tempo_espera;
    m->somatorio_tempo_espera_ao_quadrado += tempo_espera*tempo_espera;
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
        c = fila->array[fila->fim];
        fila->fim = fila_posicao_anterior(fila->fim, fila);
    }
    else {
        assert(0 && "disciplina invalida");
    }
    return c;
}

void
avanca_simulacao(Fila *fila, ListaEventos *eventos)
{
    for(int i = 0; i < eventos->tamanho_lista; ++i) {
        Evento e = eventos->array[i];
        atualiza_metricas_clientes_na_espera(e.tempo, fila->Nq, &fila->metricas);
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
            atualiza_metricas_tempo_espera(cs, &fila->metricas);
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

void
reseta_lista_de_eventos(ListaEventos *lista)
{
    bzero(lista->array, sizeof(Evento)*lista->tamanho_lista);
    lista->tamanho_lista = 0;
}

int
main()
{
    ListaEventos *eventos = (ListaEventos*)calloc(1,sizeof(*eventos)); // fazendo isso pq 'eventos' vai ter q ser ponteiro depois
    Fila *fila = (Fila*)calloc(1,sizeof(*fila));

    eventos->tamanho_array = sizeof(eventos->array)/sizeof(*eventos->array);
    fila->tamanho_array = sizeof(fila->array)/sizeof(*fila->array);
    fila->disciplina = FILA_FCFS;

    // gera primeira chegada
    {
        Evento chegada = {
            .tempo = 0,
            .tipo = EVENTO_CHEGADA_NA_FILA,
        };
        insere(chegada, eventos);
    }
    for(double t = 0; fila->metricas.coletas_clientes < 3500*1000;) {
        Evento servico = {
            .tempo = t + gera_amostra(DISTRIBUICAO_EXPONENCIAL, 1),
            .tipo = EVENTO_SERVICO_COMPLETO,
        };

        Evento chegada = {0};
        for(;;) {
            chegada.tempo = t + gera_amostra(DISTRIBUICAO_EXPONENCIAL, 0.001),
            chegada.tipo = EVENTO_CHEGADA_NA_FILA,
            t = chegada.tempo;
            if(servico.tempo < t) {
                break;
            }
            insere(chegada, eventos);
        }
        insere(servico, eventos);
        avanca_simulacao(fila, eventos);
        reseta_lista_de_eventos(eventos);
        insere(chegada, eventos);
    }
}
