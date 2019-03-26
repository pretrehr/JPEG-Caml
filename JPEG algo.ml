let round x =
	int_of_float (floor(x +. 0.5));;

let float x =
   float_of_int x;;

let pi = 4.*.atan(1.);;

let f_bloc fonction mat_bloc = (*On applique "fonction" à chaque bloc de mat_bloc*)
	let h = vect_length mat_bloc in
	let l = vect_length mat_bloc.(0) in
	let out = make_matrix h l (fonction (mat_bloc.(0).(0))) in
	for i = 0 to (h-1) do
		for j = 0 to (l-1) do
			out.(i).(j) <- fonction (mat_bloc.(i).(j))
		done;
	done;
	out;;

let produit_matrice a b =
	let m = vect_length b.(0) and n = vect_length a in
	let p = make_matrix n m 0. in
	let somme = ref 0. in
	for i = 0 to n-1 do
		for j = 0 to m-1 do
			for k = 0 to (vect_length b)-1 do
				somme := !somme +. a.(i).(k)*.b.(k).(j);
			done;
			p.(i).(j) <- !somme;
			somme := 0.;
		done;
	done;
	p;;

let transposee a =
	let n = vect_length a and m = vect_length a.(0) in
	let p = make_matrix m n a.(0).(0) in
	for i=0 to m-1 do
		for j = 0 to n-1 do
			p.(i).(j) <- a.(j).(i)
		done;
	done;
	p;;

let C x = match x with
	|0 -> 1./.sqrt(2.)
	|x when x>0 -> 1.
	|_ -> failwith "C";;

let matc n =  (*Matrice cosinus utile dans le calcul de la DCT*)
	let matrice_c = make_matrix n n (1./.sqrt(float n)) in
	for i = 1 to n-1 do
		for j = 0  to n-1 do
			matrice_c.(i).(j) <- sqrt(2./.(float(n)))*.cos((float((2*j+1)*i)*.pi)/.(float(2*n)));
		done;
	done;
	matrice_c;;

let dct_rapide matrice =
   let m = matc 8 in
   let mat = produit_matrice (produit_matrice m matrice) (transposee m) in
   let out = make_matrix 8 8 0 in
   for i = 0 to 7 do
      for j = 0 to 7 do
         out.(i).(j) <- round (mat.(i).(j))
      done;
   done;
   out;;

let dct_inverse_rapide matrice =
   let m = matc 8 in
   let out = make_matrix 8 8 0. in
   for i = 0 to 7 do
      for j = 0 to 7 do
         out.(i).(j) <- float (matrice.(i).(j))
      done;
   done;
   produit_matrice (produit_matrice (transposee m) out) m;;

let matrice_de_quantification qualite =
	let out = make_matrix 8 8 0 in
	for i = 0 to 7 do
		for j = 0 to 7 do
			out.(i).(j) <- 1 + (i+j+1)*qualite
		done;
	done;
	out;;

let matrice_quantifiee matrice qualite =
	let out = make_matrix 8 8 0 in
	let matrice_quantification = matrice_de_quantification qualite in
	for i = 0 to 7 do
		for j = 0 to 7 do
			out.(i).(j) <- (matrice.(i).(j))/(matrice_quantification.(i).(j))
		done;
	done;
	out;;

let matrice_bloc_quantifiee mat_bloc qualite =
	let h = vect_length mat_bloc in
	let l = vect_length mat_bloc.(0) in
	let out = make_matrix h l (make_matrix 8 8 0) in
	for i = 0 to (h-1) do
		for j = 0 to (l-1) do
			out.(i).(j) <- matrice_quantifiee (mat_bloc.(i).(j)) qualite
		done;
	done;
	out;;

let matrice_dequantifiee matrice qualite =
	let out = make_matrix 8 8 0 in
	let matrice_quantification = matrice_de_quantification qualite in
	for i = 0 to 7 do
		for j = 0 to 7 do
			out.(i).(j) <- (matrice.(i).(j))*(matrice_quantification.(i).(j))
		done;
	done;
	out;;

let matrice_bloc_dequantifiee mat_bloc qualite = 
	let h = vect_length mat_bloc in
	let l = vect_length mat_bloc.(0) in
	let out = make_matrix h l (make_matrix 8 8 0) in
	for i = 0 to (h-1) do
		for j = 0 to (l-1) do
			out.(i).(j) <- matrice_dequantifiee (mat_bloc.(i).(j)) qualite
		done;
	done;
	out;;

let erreur image qualite =
   let out = make_matrix 8 8 0. in
   let apres = dct_inverse_rapide (matrice_dequantifiee (matrice_quantifiee (dct_rapide image) qualite) qualite) in
   for i = 0 to 7 do
      for j = 0 to 7 do
         out.(i).(j) <- image.(i).(j) -. apres.(i).(j)
      done;
   done;
   out;;

let linearisation matrice =
	let tab = make_vect 64 matrice.(0).(0) in
	let x = ref 0 and y = ref 0 in
	for i=1 to 7 do
		if (i mod 2)=0 then 
			(x := i; y:= 0)
		else (x:= 0; y:= i);
		for k = (i*(i+1)/2) to ((i+2)*(i+1)/2-1) do
			tab.(k) <- matrice.(!x).(!y);
			if (i mod 2)=0 then
				(decr x;	incr y)
			else (decr y;	incr x)
		done;
	done;
	for i=1 to 7 do
		if (i mod 2)=0 then (x := 7-i; y:= 7) 
		else (x := 7; y := 7-i);
		for k = (63-(i*(i+1)/2)) downto (63-((i+2)*(i+1)/2-1)) do
			tab.(k) <- matrice.(!x).(!y);
			if (i mod 2)=0 then
				(incr x; decr y)
			else (decr x; incr y)
		done;
	done;
	tab.(63) <- matrice.(7).(7);
	tab;;

let delinearisation tab =
	let out = make_matrix 8 8 tab.(0) in
	let x = ref 0 and y = ref 0 in
	for i = 1 to 7 do
		if (i mod 2)=0 then 
			(x := i; y:= 0)
		else (x:= 0; y:= i);
		for k = (i*(i+1)/2) to ((i+2)*(i+1)/2-1) do
			out.(!x).(!y) <- tab.(k) ;
			if (i mod 2)=0 then
				(decr x;	incr y)
			else (decr y; incr x)
		done;
	done;
	for i=1 to 7 do
		if (i mod 2)=0 then 
			(x := 7-i; y:= 7) 
		else (x := 7; y := 7-i);
		for k = (63-(i*(i+1)/2)) downto (63-((i+2)*(i+1)/2-1)) do
			out.(!x).(!y) <- tab.(k);
			if (i mod 2)=0 then
				(incr x; decr y)
			else (decr x; incr y)
		done;
	done;
	out.(7).(7) <- tab.(63);
	out;;

let reduction tab = (*Codage RLE d'un tableau ([|5;3;0;0;0;0;2|]->[|5;3;0;4;2|]) *)
	let l = list_of_vect tab in
	let rec aux liste compteur0 = match liste with (*liste est précédée de "compteur0" zéros*)
		|[] when compteur0 <> 0 -> compteur0::[]
		|[] -> []
		|x::queue when x=0 -> if compteur0=0 then 
										 0::(aux queue (compteur0+1))
									 else (aux queue (compteur0+1))
		|x::queue -> if compteur0 <> 0 then 
							 compteur0::x::(aux queue 0)
						 else x::(aux queue 0)
	in vect_of_list (aux l 0);;

let rec liste_de_0 n = match n with (*génère une liste de n zéros*)
	|0->[]
	|n-> 0::(liste_de_0 (n-1));;

let dereduction tab = (*effectue l'opération inverse du codage RLE*)
	let l = list_of_vect tab in
	let rec aux l = match l with
	|[] -> []
	|0::x::queue -> (liste_de_0 x)@(aux queue)
	|x::queue -> x::(aux queue)
	in vect_of_list (aux l);;

type ' a arbre =
    | Feuille of ' a
    | Noeud of ' a arbre * ' a arbre;;

let rec occ caractere liste = match liste with (*ajoute une occurence dans la liste d'occurences pour le caractère considéré*)
		|[] -> [1,caractere]
		|(nx,x)::queue when x=caractere -> (nx+1,caractere)::queue
		|(nx,x)::queue -> (nx,x)::(occ caractere queue);;

let tri liste = 
	sort__sort prefix <= liste;;

let occurences tab = (*renvoie la liste des occurences de chaque élément de tab*)
	let n = vect_length tab in
	let list = ref [] in
	for i = 0 to (n-1) do
		list := occ tab.(i) !list;
	done;
	!list;;

let rec insere x liste = match liste with
    | [] -> [x]
    | y :: queue -> if x < y then x :: y :: queue
        else y :: (insere x queue);;

let rec fusionne liste = match liste with
    | [] -> failwith "fusionne"
    | [nx, x] -> x
    | (nx, x) :: (ny, y) :: queue -> fusionne (insere (nx + ny, Noeud (y, x)) queue);;

let huffman liste =
    let l = map (function (n, x) -> (n, Feuille x)) (tri liste) in
    fusionne l;;

let arbre_de_huffman tab =
    huffman (occurences tab);;

let retourne liste =
  let rec aux acc l =
    if l = []
    then acc
    else aux (hd(l) :: acc) (tl l)
  in aux [] liste;;

let table_equivalence arbre = (*renvoie l'ensemble des caractères de l'arbre et leur chemin d'accès dans l'arbre i.e. la table binaire*)
	let rec aux chemin arbre = match arbre with
		|Feuille x -> [x,retourne chemin]
		|Noeud(x,y) -> (aux (false::chemin) x)@(aux (true::chemin) y)
	in aux [] arbre;;

let code vect = (*renvoie l'arbre de Huffman et la liste de booléen associé à vect*)
    let a = arbre_de_huffman vect in
    let codes_car = table_equivalence a in
    let n = vect_length vect in
    let rec aux liste compteur = match compteur with
        | - 1 -> (a, liste)
        | _ -> aux ((assoc vect.(compteur) codes_car) @ liste) (compteur - 1)
    in aux [] (n - 1);;

let rec lire a w = match (a, w) with (*renvoie la première feuille atteinte par le mot w dans l'arbre a, ainsi que le reste du mot qui n'a pas encore été lu*)
    | (Feuille x, l) -> (x, l)
    | (Noeud (x, y), []) -> failwith "lire"
    | (Noeud (x, y), d :: queue) -> if d then lire y queue
        else lire x queue;;

let decodage (a,w) = (*renvoie le vecteur représenté par le mot binaire w dans l'arbre a*)
	let rec decode arbre mot_bin = match mot_bin with
   	 | [] -> []
	    | m -> let (x, m1) = lire a m in
   	     x::(decode arbre m1)
	in vect_of_list(decode a w);;

let list_of_string chaine =
	let n = string_length chaine in 
	let tab = make_vect n 0 in
	for i = 0 to n-1 do
		tab.(i) <- int_of_char chaine.[i]
	done;
	let liste = list_of_vect tab in
	let list = map (function 49 -> true | _ -> false) liste in
	list;;

let taille_bit arbre = (*renvoie le nb de bits qu'occupe un arbre de Huffman*)
	let table = table_equivalence arbre in
	let rec aux table = match table with
		|[] -> 0
		|(x,chemin)::queue -> 2*8+(list_length chemin) + (aux queue)
	in aux table;;

let taille_image destination = (*renvoie les dimensions d'une image située dans le système*)
	let canal = open_in destination in (*open_in ouvre en lecture le fichier dont le nom est passé en paramètre*)
	let largeur = ref 0 in
	seek_in canal 18;
	for i = 0 to 3 do
		largeur := !largeur + (input_byte canal)*int_of_float(16.**float(i*2))
	done;
	let hauteur = ref 0 in
	for j = 0 to 3 do
		hauteur := !hauteur + (input_byte canal)*int_of_float(16.**float(j*2))
	done;
	close_in canal;
	!largeur,!hauteur;;

#open "graphics";;

let matrice_de_fichier destination = (*renvoie la matrice des couleurs de l'image*)
    let canal=open_in_bin destination in
    let (largeur,hauteur)=taille_image destination in    
    let (m:color vect vect)=make_matrix hauteur largeur 0 in    
    let b=ref 0 and v=ref 0 and r=ref 0 and col=ref 0 in    
    seek_in canal 54;
    for i=(hauteur-1) downto 0 do
        for j=0 to (largeur-1) do            
        		b:=input_byte canal;            
        		v:=input_byte canal;            
        		r:=input_byte canal;            
        		col:=(rgb (!r) (!v) (!b));            
        		m.(i).(j)<-(!col)        
        done
    done;    
    close_in canal;    
    m;;

let redim destination = (*On redimensionne la matrice pour rendre l et h multiples de 8*)
   let matrice = matrice_de_fichier destination in
   let (l, h) = taille_image destination in
   let resteh = (h mod 8) in
   let restel = (l mod 8) in
   if resteh = 0 && restel = 0 then matrice
   else let out = make_matrix (h + 8 - resteh) (l + 8 - restel) white in
      for i = 0 to (h + 8 - resteh - 1) do
         for j = 0 to (l + 8 - restel - 1) do
            if (i < h && j < l) then out.(i).(j) <- matrice.(i).(j)
         done;
      done;
      out;;

let colors color = (*renvoie le triplet [|r;g;b|] à partir d'un coefficient de type color*)
	let b = color mod 256 in
	let g = ((color-b) mod (256*256))/256 in
	let r = ((color-b-g) mod (256*256*256))/(256*256) in
	[|r;g;b|];;

let matrices_rgb image = (* renvoie les matrices r, g et b à partir de la matrice de type color*)
	let hauteur = vect_length image in
	let longueur = vect_length image.(0) in
	let matr = make_matrix hauteur longueur image.(0).(0)
	and matg = make_matrix hauteur longueur image.(0).(0)
	and matb = make_matrix hauteur longueur image.(0).(0) in
	let couleurs = ref (colors image.(0).(0)) in
	for i=0 to (hauteur-1) do
		for j = 0 to (longueur -1) do
			couleurs := (colors image.(i).(j));
			matr.(i).(j) <- !couleurs.(0);
			matg.(i).(j) <- !couleurs.(1);
			matb.(i).(j) <- !couleurs.(2);
		done;
	done;
	[|matr;matg;matb|];;

let matrices_ycbcr image = (*convertit les matrices R,G et B en Y, Cb et Cr*)
	let hauteur = vect_length image in
	let longueur = vect_length image.(0) in
	let tab_rgb = matrices_rgb image in
	let maty = make_matrix hauteur longueur 0.
	and matcb = make_matrix hauteur longueur 0.
	and matcr = make_matrix hauteur longueur 0. in
	for i = 0 to (hauteur-1) do
		for j = 0 to (longueur-1) do
			let r = ref tab_rgb.(0).(i).(j)
			and g = ref tab_rgb.(1).(i).(j)
			and b = ref tab_rgb.(2).(i).(j) in
			maty.(i).(j) <- 0.299*.float(!r)+.0.587*.float(!g)+.0.114*.float(!b);
			matcb.(i).(j) <- -0.1687*.float(!r)-.0.3313*.float(!g)+.0.5*.float(!b)+.128.;
			matcr.(i).(j) <- 0.5*.float(!r)-.0.4187*.float(!g)-.0.0813*.float(!b)+.128.
		done;
	done;
	[|maty;matcb;matcr|];;

let new_matrices_rgb maty matcb matcr = (*convertit les matrices Y, Cb et Cr en R, G et B*)
	let hauteur = vect_length maty in
	let longueur = vect_length maty.(0) in
	let matr = make_matrix hauteur longueur 0 
	and matg = make_matrix hauteur longueur 0
	and matb = make_matrix hauteur longueur 0 in
	for i = 0 to (hauteur-1) do
		for j = 0 to (longueur-1) do
			let y = ref maty.(i).(j)
			and cb = ref matcb.(i).(j)
			and cr = ref matcr.(i).(j) in
			matr.(i).(j) <- max (min (round (!y+.1.402*.(float(!cr)-.128.))) 255) 0;
			matg.(i).(j) <- max (min (round (!y-.0.34414*.(float(!cb)-.128.)-.0.71414*.(float(!cr)-.128.))) 255) 0;
			matb.(i).(j) <- max (min (round (!y+.1.772*.(float(!cb)-.128.))) 255) 0
		done;
	done;
	[|matr;matg;matb|];;

let regroupement_rgb matr matg matb = (*renvoie la matrice de type color associée aux matrices R, G et B*)
	let hauteur = vect_length matr in
	let longueur = vect_length matr.(0) in
	let (out: color vect vect) = make_matrix hauteur longueur 0 in
	for i = 0 to (hauteur-1) do
		for j = 0 to (longueur-1) do
			out.(i).(j) <- rgb (matr.(i).(j)) (matg.(i).(j)) (matb.(i).(j))
		done;
	done;
	out;;

let share_line matrice_ligne cote = (* partage une matrice de hauteur cote en plusieurs matrices de largeur cote *)
   let rec partager matrice compteur =
      if compteur = ((vect_length matrice.(0)) / cote) then []
      else let mat = ref (make_matrix cote cote matrice.(0).(0)) in
         for j = 0 to (cote - 1) do
            for k = 0 to (cote - 1) do
               !mat.(j).(k) <- matrice.(j).(cote * compteur + k)
            done;
         done;
         !mat :: (partager matrice (compteur + 1))
   in vect_of_list (partager matrice_ligne 0);;

let share matrice cote = (*partage une matrice en matrices de taille cote*)
   if cote > vect_length matrice || cote > vect_length matrice.(0) then failwith "partage impossible";
   let rec aux_share mat compteur =
      if (compteur = ((vect_length mat) / cote)) then []
      else (share_line (sub_vect mat (cote * compteur) cote) cote) :: (aux_share mat (compteur + 1))
   in vect_of_list (aux_share matrice 0);;

let regroupement_blocs mat_bloc = (*opération inverse de la fonction share*)
	let h = vect_length mat_bloc in
	let l = vect_length mat_bloc.(0) in
	let cote = vect_length mat_bloc.(0).(0) in
	let out = make_matrix (h*cote) (l*cote) mat_bloc.(0).(0).(0).(0) in
	for i = 0 to (h-1) do
		for j = 0 to (l-1) do
			for y = 0 to (cote-1) do
				for x = 0 to (cote-1) do
					out.(cote*i+y).(cote*j+x) <- mat_bloc.(i).(j).(y).(x)
				done;
			done;
		done;
	done;
	out;;

let moyenne bloc =
	round ((bloc.(0).(0)+.bloc.(1).(0)+.bloc.(0).(1)+.bloc.(1).(1))/.4.);;

let matrice_reduite matrice = (*réduit la matrice en effectuant la moyenne de chaque bloc 2x2*)
	let tab = share matrice 2 in
	f_bloc moyenne tab;;

let agrandissement matrice = (*opération inverse de la fonction matrice_reduite*)
	let h = vect_length matrice in
	let l = vect_length matrice.(0) in
	let out = make_matrix (h*2) (l*2) 0(*.*) in
	for x = 0 to (h-1) do
		for y = 0 to (l-1) do
			let valeur = ref matrice.(x).(y) in
			out.(2*x).(2*y) <- !valeur;
			out.(2*x+1).(2*y) <- !valeur;
			out.(2*x).(2*y+1) <- !valeur;
			out.(2*x+1).(2*y+1) <- !valeur;
		done;
	done;
	out;;

matrice_reduite [|[|1.;2.;3.;4.|];[|5.;6.;7.;8.|];[|9.;10.;11.;12.|];[|13.;14.;15.;16.|]|];;



let jpeg chemin niveau =
   (*etapes de la compression :
	-partage de Y en blocs de 8x8
	-transformée en cosinus discrète
	-quantification
	-linearisation
	-RLE sur les 0
	-codage de Huffman
	-moyennage des matrices Cb,Cr*)
   let a = sys__time () in
   let largeur, hauteur = taille_image chemin in
   let y = (matrices_ycbcr (redim chemin)).(0) in
   let cb = (matrices_ycbcr (redim chemin)).(1) in
   let cr = (matrices_ycbcr (redim chemin)).(2) in
   let out = f_bloc code (f_bloc reduction (f_bloc linearisation (matrice_bloc_quantifiee (f_bloc dct_rapide (share y 8)) niveau))), (matrice_reduite cb, matrice_reduite cr) in
   print_float (sys__time () -. a); print_newline ();
   out;;

let decompression (yrle, (cb, cr)) niveau =
   let a = sys__time () in
   let ybloc = f_bloc delinearisation (f_bloc dereduction (f_bloc decodage yrle)) in
   let yapres = f_bloc dct_inverse_rapide (matrice_bloc_dequantifiee ybloc niveau) in
   let cbbis = agrandissement cb in
   let crbis = agrandissement cr in
   let ybis = regroupement_blocs yapres in
   let new_mat = (regroupement_rgb (new_matrices_rgb ybis cbbis crbis).(0) (new_matrices_rgb ybis cbbis crbis).(1) (new_matrices_rgb ybis cbbis crbis).(2)) in
   let hauteur = vect_length ybis in
   let largeur = vect_length ybis.(0) in
   let chaine = (string_of_int (min 1350 largeur)) ^ "x" ^ (string_of_int (min 690 hauteur)) in
   open_graph chaine;
   print_float (sys__time () -. a);
   draw_image (make_image new_mat) 0 0;;

let taille_bit arbre =
   let table = table_equivalence arbre in
   let rec aux table = match table with
      | [] -> 0
      | (x, chemin) :: queue -> 2 * 8 + (list_length chemin) + (aux queue)
   in aux table;;

let taille_bit_y mat = 
	let h = vect_length mat and l = vect_length mat in
	let out = ref 0 in
	for i = 0 to (h-1) do
		for j = 0 to (l-1) do
			out:=!out +(taille_bit(fst mat.(i).(j)))+ (list_length (snd (mat.(i).(j))))
		done;
	done;
	!out;;

let taille_bit_cbcr mat =
	let h = vect_length mat and l = vect_length mat in
	8*h*l;;

let taille_compressee (yrle, (cb, cr)) =
	print_newline();
	print_int ((taille_bit_y yrle + 2*(taille_bit_cbcr cb))/8000); print_string " ko";;

redim "F:/MPE/TIPE/Bitmap pdf/pdf7.bmp";;
taille_image "F:/MPE/TIPE/Bitmap pdf/pdf7.bmp";;
matrice_de_fichier "F:/MPE/TIPE/Bitmap pdf/pdf8.bmp";;
jpeg2 "F:/MPE/TIPE/Bitmap pdf/pdf.bmp" 3;;
decompression (jpeg "F:/MPE/TIPE/tiger.bmp" 0) 0;;
matrice_reduite;;
linearisation;;
matrice_bloc_quantifiee;;
f_bloc linearisation (matrice_bloc_quantifiee (f_bloc dct_rapide (share ((matrices_ycbcr (matrice_de_fichier "F:/MPE/TIPE/Bitmap pdf/pdf8.bmp")).(0)) 8)) 2);;
matrices_rgb (matrice_de_fichier "F:/MPE/TIPE/Bitmap pdf/pdf8.bmp");;
taille_image "F:/MPE/TIPE/Bitmap pdf/pdf8.bmp";;
matrices_ycbcr;;
share;;
f_bloc reduction;;


decompression (jpeg2 "F:/MPE/TIPE/Hiver.bmp" 3) 3;;


taille_compressee (jpeg "F:/MPE/TIPE/tiger.bmp" 1000);;
