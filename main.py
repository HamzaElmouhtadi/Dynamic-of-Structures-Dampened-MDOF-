import numpy as np
from tkinter import Tk, Label, StringVar, Button, Entry
from scipy.linalg import eig
from calculator import Calc

window = Tk()
screen_width = window.winfo_screenwidth()
screen_height = window.winfo_screenheight()
Size = (screen_height,screen_width)
window.title("Dynamique des Structures")
window.geometry("800x500+120+120")
window.resizable(True, True)
window.state("zoomed")
N = 0


text_var_k = []
text_var_b = []
text_var_m = []
entries_k = []
entries_b = []
entries_m = []
text_var_f = []
entries_f = []
def update_K(text) :
     Lbl_K.config(text="K\n"+text)
def update_M(text) :
     Lbl_M.config(text="M\n"+text)
def update_S(text) :
     Lbl_S.config(text="Ʌ\n"+text)
def update_X(text) :
     Lbl_X.config(text="X\n"+text)
def update_B(text) :
     Lbl_B.config(text="B\n"+text)
def update_Freq(text) :
     Lbl_F.config(text="F Fondamentale (Hz)\n"+text)
def Test():
     N=3
     M = np.array([[6.0000, 0, 0], [0, 0.5000, 0], [0, 0, 0.0800]])
     K = np.array([[1013000, -13000, 0], [-13000, 44000, -11000], [0, -11000, 11500]])
     B = np.array([[0.4, -0.4, 0], [-0.4, 1.1, 0], [0, 0, 0.7]])
     F = np.array([0, 0, 1000])
# Appel de la fonction Calc avec les valeurs assignées
     Calc(3, M, K, F, B)


      
def get_mat():
     N=int(Entry_N.get())
     m=np.array(np.zeros((N),dtype=float))  #Les masses 
     k=np.array(np.zeros((N,N),dtype=float))  #Les raideurs
     M=np.array(np.zeros((N,N),dtype=float))  #Matrise Masse 
     K=np.array(np.zeros((N,N),dtype=float))  #Matrice raideur
     b = np.array(np.zeros((N, N), dtype=float))  # Matrice Ammortissement
     B = np.array(np.zeros((N, N), dtype=float))  # Matrice Ammortissement
     S=np.array(np.zeros((N,N),dtype=float))  #Matrice Spectrale
     X=np.array(np.zeros((N,N),dtype=float))  #Matrice Modale
     F=np.array(np.zeros((N,1),dtype=float))  #Vecteur Force
     temps = 0
     for i in range(0,N) :
          m[i] = float(text_var_m[i][0].get())
          F[i] = float(text_var_f[i][0].get())
          for j in range(0,i+1) :
               val = float(text_var_k[i][j].get())
               val2 = float(text_var_b[i][j].get())
               k[i][j] = val
               k[j][i] = val
               b[i][j] = val2
               b[j][i] = val2

     M = np.diag(m)
     K = np.diag(np.sum(k, axis=0)) - np.triu(k, 1) - np.tril(k, -1)
     B = np.diag(np.sum(b, axis=0)) - np.triu(b, 1) - np.tril(b, -1)

     #resolution de K*x = lambda*M*x

     valeurs_propres,X = eig(K, M)
     S=np.diag(np.round(valeurs_propres,3))

     min_valeur_propre_index = np.argmin(valeurs_propres)
     lambda_min = valeurs_propres[min_valeur_propre_index]
     update_K(str(np.matrix(K)))
     update_M(str(np.matrix(M)))
     update_S(str(np.matrix(S)))
     update_X(str(np.matrix(X)))
     update_B(str(np.matrix(B)))
     update_Freq(str(round(np.sqrt(lambda_min)/(2*np.pi),6)))
     Calc(N,M,K,F,B)



Entry_N = Entry(window,textvariable = N,width=10)
Entry_N.place(x = 190,y = 0)
Lbl_K = Label(window, text="K\n", font=('slant', 20, 'bold'), )
Lbl_K.place(x=250, y=300)
Lbl_M = Label(window, text="M\n", font=('Helvetica', 20, 'bold'),  )
Lbl_M.place(x=0, y=300)
Lbl_F = Label(window, text="F Fondamentale (Hz)\n", font=('arial', 20, 'bold'),)
Lbl_F.place(x=500, y=300)
Lbl_S = Label(window, text="Ʌ\n", font=('Helvetica', 20, 'bold'),  )
Lbl_S.place(x=800, y=300)
Lbl_B = Label(window, text="B\n", font=('Helvetica', 20, 'bold'),  )
Lbl_B.place(x=0, y=500)
Lbl_X = Label(window, text="X\n", font=('arial', 20, 'bold'),)
Lbl_X.place(x=500, y=500)
Label(window, text="Entrer la taille de la matrice  :", font=('slant', 8, 'bold'),).place(x=20, y=0)
Label(window, text="Entrer les raideurs:", font=('slant', 10, 'bold'),).place(x=20, y=20)
Label(window, text="M:", font=('slant', 10, 'bold'),).place(x=250, y=20)
Label(window, text="F", font=('slant', 10, 'bold'),).place(x=350, y=20)
Label(window, text="B", font=('slant', 10, 'bold'),).place(x=450, y=20)
def draw_mat():
     x2 = 0
     y2 = 0
     N=int(Entry_N.get())
     for i in range(N): 
          text_var_k.append([])
          text_var_b.append([])
          text_var_m.append([])
          entries_k.append([])
          entries_b.append([])
          entries_m.append([])
          text_var_f.append([])
          entries_f.append([])
          for j in range(i+1):
               text_var_k[i].append(StringVar())
               text_var_b[i].append(StringVar())
               text_var_m[i].append(StringVar())
               text_var_f[i].append(StringVar())
               entries_k[i].append(Entry(window, textvariable=text_var_k[i][j],width=3))
               entries_k[i][j].place(x=60 + x2, y=50 + y2)
               entries_b[i].append(Entry(window, textvariable=text_var_b[i][j], width=3))
               entries_b[i][j].place(x=450 + x2, y=50 + y2)
               entries_m[i].append(Entry(window, textvariable=text_var_m[i][0],width=3))
               entries_m[i][0].place(x=250, y=50 + y2)
               entries_f[i].append(Entry(window, textvariable=text_var_f[i][0], width=3))
               entries_f[i][0].place(x=350, y=50 + y2)
               x2 += 30
          y2 += 30
          x2 = 0

def reset() :
     update_M("")
     update_K("")
     update_B("")
     update_Freq("")
     update_S("")
     update_X("")
     try:
          for i in range(len(entries_m)):
               entries_m[i][0].place_forget()
               entries_f[i][0].place_forget()
               for j in range(i+1):
                    entries_k[i][j].place_forget()
                    entries_b[i][j].place_forget()
     except IndexError:
          pass 
     entries_k.clear()
     entries_b.clear()
     entries_m.clear()
     text_var_k.clear()
     text_var_b.clear()
     text_var_m.clear()
     entries_f.clear()
     text_var_f.clear()



button= Button(window,text="Calculate", bg='bisque3', width=15, command=get_mat)
button.place(x=200,y=200)
button= Button(window,text="Test", bg='bisque3', width=15, command=Test)
button.place(x=300,y=200)
button2= Button(window,text="Reset", bg='bisque3', width=15, command=reset)
button2.place(x=350,y=0)
button3= Button(window,text="Build", bg='bisque3', width=15, command=draw_mat)
button3.place(x=230,y=0)



window.mainloop()