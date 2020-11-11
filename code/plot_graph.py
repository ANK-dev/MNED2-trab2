###############################################################################
#         Métodos Numéricos para Equações Diferenciais II -- Trabalho 2       #
#                           Ariel Nogueira Kovaljski                          #
###############################################################################

import matplotlib.pyplot as plt


x0, y0 = [], []
x1, y1 = [], []
x2, y2 = [], []
x3, y3 = [], []
methods = {0: 'FTBS', 1: 'Lax-Friedrichs', 2: 'Lax-Wendroff', 3: 'Beam-Warming'}
parameters = {}

def main():
    # Abre arquivos para leitura
    # Arquivo FTBS
    try:
        f0 = open('./results/results0.txt', 'r')
    except OSError as err:
        print("Erro ao abrir o arquivo 'results0.txt':", err)

    # Arquivo Lax-Friedrichs
    try:
        f1 = open('./results/results1.txt', 'r')
    except OSError as err:
        print("Erro ao abrir o arquivo 'results1.txt':", err)

    # Arquivo Lax-Wendroff
    try:
        f2 = open('./results/results2.txt', 'r')
    except OSError as err:
        print("Erro ao abrir o arquivo 'results2.txt':", err)

    # Arquivo Beam-Warming
    try:
        f3 = open('./results/results3.txt', 'r')
    except OSError as err:
        print("Erro ao abrir o arquivo 'results3.txt':", err)

    # Extrai informações dos arquivos
    file_format(f0, x0, y0)
    file_format(f1, x1, y1)
    file_format(f2, x2, y2)
    file_format(f3, x3, y3)

    # Plota os gráficos de cada método
    plot_graph(x0, y0, 0)
    plot_graph(x1, y1, 1)
    plot_graph(x2, y2, 2)
    plot_graph(x3, y3, 3)

    # Exibe os gráficos
    plt.show()


def file_format(f, x, y):
    for line_number, line in enumerate(f):
        if line_number < 8:
            # Adiciona parâmetros da simulação em um dicionário
            parameters[line.split('=')[0]] = (line.split('=')[1]).split('\n')[0]
        elif line_number == 8:
            # Pula linha separadora
            pass
        else:
            # Separa valores na lista `x` e na lista `y`
            x.append( int(line.split(',')[0]) )
            y.append( float( (line.split(',')[1]).split('\n')[0] ) )


def plot_graph(x, y, id):
    # Configura e exibe o gráfico
    colors = {0: 'red', 1: 'green', 2: 'blue', 3: 'yellow'}
    fig,ax = plt.subplots()
    fig.set_size_inches(8, 7)   # Size of the window (1in = 100px)
    ax.grid(True)

    plt.suptitle(f"Método {methods[id]}: Concentração X Célula")
    plt.title(rf"$nx = {parameters['nx']}$, " 
              rf"$\Delta t = {parameters['Delta_t']}$, "
              rf"$\Delta x = {parameters['Delta_x']}$, "
              rf"$t_{{final}} = {parameters['t_final']}$, "
              rf"$\bar{{u}} = {parameters['u_bar']}$, "
              rf"$\alpha = {parameters['alpha']}$, "
              rf"$c_{{ini}} = {parameters['c_ini']}$, "
              rf"$c_{{inj}} = {parameters['c_inj']}$", fontsize=8)

    plt.xlabel("índice célula (i)")
    plt.ylabel("concentração (Q)")
    plt.plot(x,y,'ko-', markerfacecolor=colors[id], markeredgecolor='k')

if __name__ == "__main__":
    main()
