import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from scipy.signal import welch
import scipy.io

class Application(tk.Tk):
    def __init__(self):
        super().__init__()

        # Initialisation de la fenêtre principale
        self.title("Analyse de Signal")
        self.geometry("1200x800")

        # Variables globales pour stocker les signaux
        self.signal_original = None
        self.signal_actuel = None

        # Création des widgets de l'interface
        self.create_widgets()

    def create_widgets(self):
        # Frame principale pour le graphique
        self.frame_graph = ttk.Frame(self)
        self.frame_graph.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Initialisation de la figure matplotlib
        self.figure = Figure(figsize=(5, 4), dpi=100)
        self.ax = self.figure.add_subplot(111)

        # Canvas pour intégrer matplotlib avec tkinter
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.frame_graph)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Frame pour les boutons
        self.frame_buttons = ttk.Frame(self)
        self.frame_buttons.pack(side=tk.RIGHT, fill=tk.Y)

        # Boutons et zones de saisie
        self.create_button("Charger signal", self.load_signal)
        self.create_label_button_pair("RSB (dB):", "Bruiter signal", self.noise_signal)
        self.create_button("Réinitialiser signal", self.reload_signal)
        self.create_button("Afficher périodogramme", self.show_periodogram)
        self.create_button("Afficher profil", self.show_profile)
        self.create_label_button_pair("N:", "Afficher découpage", self.show_cut)
        self.create_button("Afficher résidus", self.show_residue)
        self.create_button("Afficher F²(N)", self.show_F2N)

    def create_button(self, text, command):
        button = ttk.Button(self.frame_buttons, text=text, command=command)
        button.pack(pady=10, padx=5, fill=tk.X)

    def create_label_button_pair(self, label_text, button_text, command):
        frame = ttk.Frame(self.frame_buttons)
        frame.pack(pady=10, padx=5, fill=tk.X)

        label = ttk.Label(frame, text=label_text)
        label.pack(side=tk.LEFT, padx=5)

        entry = ttk.Entry(frame)
        entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        button = ttk.Button(frame, text=button_text, command=lambda: command(entry))
        button.pack(side=tk.LEFT, padx=5)

    def load_signal(self):
        file_path = filedialog.askopenfilename(filetypes=[("MAT files", "*.mat")])
        if file_path:
            try:
                # Chargement du fichier .mat
                mat = scipy.io.loadmat(file_path)
                
                # Trouver les variables qui ne commencent pas par '__'
                variable_names = [name for name in mat.keys() if not name.startswith('__')]
                
                if 'data' in variable_names:  # Pour le fichier data_Weierstrass.mat
                    # On prend la première régularité (indice 0)
                    self.signal_original = mat['data'][0,0].flatten()
                else:
                    # Essayer de trouver une variable qui pourrait être le signal
                    for var_name in variable_names:
                        if isinstance(mat[var_name], np.ndarray):
                            self.signal_original = mat[var_name].flatten()
                            break
                    else:
                        messagebox.showerror("Erreur", f"Pas de signal trouvé dans les variables : {variable_names}")
                        return

                # Copier le signal original comme signal actuel
                self.signal_actuel = self.signal_original.copy()
                
                # Afficher le signal
                self.ax.clear()
                self.ax.plot(self.signal_actuel)
                self.ax.set_title("Signal Chargé")
                self.ax.grid(True)
                self.canvas.draw()

            except Exception as e:
                messagebox.showerror("Erreur de chargement", f"Erreur lors du chargement du fichier : {str(e)}")

    def noise_signal(self, entry):
        if self.signal_original is None:
            messagebox.showerror("Erreur", "Veuillez charger un signal d'abord.")
            return

        rsb = entry.get()
        try:
            rsb = float(rsb)
            power_signal = np.mean(self.signal_original ** 2)
            power_noise = power_signal / (10 ** (rsb / 10))
            noise = np.sqrt(power_noise) * np.random.randn(len(self.signal_original))
            self.signal_actuel = self.signal_original + noise
            self.ax.clear()
            self.ax.plot(self.signal_actuel)
            self.ax.set_title(f"Signal bruité (RSB = {rsb} dB)")
            self.canvas.draw()
        except ValueError:
            messagebox.showerror("Erreur de validation", "Veuillez entrer un RSB valide")

    def reload_signal(self):
        if self.signal_original is not None:
            self.signal_actuel = self.signal_original.copy()
            self.ax.clear()
            self.ax.plot(self.signal_actuel)
            self.ax.set_title("Signal Réinitialisé")
            self.canvas.draw()
        else:
            messagebox.showerror("Erreur", "Aucun signal à recharger")

    def show_periodogram(self):
        if self.signal_actuel is None:
            messagebox.showerror("Erreur", "Veuillez charger ou générer un signal d'abord.")
            return

        f, Pxx_den = welch(self.signal_actuel, fs=1.0, nperseg=1024)
        self.ax.clear()
        self.ax.semilogy(f, Pxx_den)
        self.ax.set_xlabel('Fréquence [Hz]')
        self.ax.set_ylabel('Densité spectrale de puissance [V**2/Hz]')
        self.ax.set_title('Périodogramme')
        self.canvas.draw()

    def show_profile(self):
        if self.signal_actuel is None:
            messagebox.showerror("Erreur", "Veuillez charger ou générer un signal d'abord.")
            return

        profil = self.cumulative_profile(self.signal_actuel)
        self.ax.clear()
        self.ax.plot(profil)
        self.ax.set_xlabel('Échantillon')
        self.ax.set_ylabel('Profil cumulatif')
        self.ax.set_title('Profil cumulatif du signal')
        self.canvas.draw()

    def show_cut(self, entry):
        if self.signal_actuel is None:
            messagebox.showerror("Erreur", "Veuillez charger ou générer un signal d'abord.")
            return

        n = entry.get()
        try:
            n = int(n)
            profil = self.cumulative_profile(self.signal_actuel)
            segments = np.array_split(profil, n)
            self.ax.clear()
            for segment in segments:
                x = np.arange(len(segment))
                coeffs = np.polyfit(x, segment, 1)
                trend = np.polyval(coeffs, x)
                self.ax.plot(x, segment)
                self.ax.plot(x, trend, linestyle='--')
            self.ax.set_xlabel('Échantillon')
            self.ax.set_ylabel('Profil cumulé')
            self.ax.set_title('Découpage et tendances locales')
            self.canvas.draw()
        except ValueError:
            messagebox.showerror("Erreur de validation", "Veuillez entrer une valeur de N valide")

    def show_residue(self):
        if self.signal_actuel is None:
            messagebox.showerror("Erreur", "Veuillez charger ou générer un signal d'abord.")
            return

        profil = self.cumulative_profile(self.signal_actuel)
        segments = np.array_split(profil, 4)
        self.ax.clear()
        for segment in segments:
            x = np.arange(len(segment))
            coeffs = np.polyfit(x, segment, 1)
            trend = np.polyval(coeffs, x)
            residue = segment - trend
            self.ax.plot(residue)
        self.ax.set_xlabel('Échantillon')
        self.ax.set_ylabel('Résidu')
        self.ax.set_title('Résidus du profil')
        self.canvas.draw()

    def show_F2N(self):
        if self.signal_actuel is None:
            messagebox.showerror("Erreur", "Veuillez charger ou générer un signal d'abord.")
            return

        profil = self.cumulative_profile(self.signal_actuel)
        n_values = np.logspace(1, 2.5, 10).astype(int)
        f2n = []

        for n in n_values:
            segments = np.array_split(profil, len(profil) // n)
            rms = np.sqrt(np.mean([np.var(segment) for segment in segments]))
            f2n.append(rms ** 2)

        self.ax.clear()
        self.ax.loglog(n_values, f2n, marker='o')
        self.ax.set_xlabel('Segments')
        self.ax.set_ylabel('F²(N)')
        self.ax.set_title('Courbe F²(N)')
        self.canvas.draw()

    def cumulative_profile(self, signal):
        signal_centre = signal - np.mean(signal)
        profil_cumulatif = np.cumsum(signal_centre)
        return profil_cumulatif

if __name__ == "__main__":
    app = Application()
    app.mainloop()
