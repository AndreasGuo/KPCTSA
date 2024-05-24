This is code for Knee Point based Improved PO. 
To run code, using
```bash
git clone https://github.com/AndreasGuo/KPIPODNA.git
cd KPIPODNA
go get
go run --fitrev
```

Other parameters may useful for comparison

| parameter | type | default | desc |
|--|--|--|--|
| vdim | Int | 20 | length of one DNA sequence |
| setsize | Int | 7 | size of target DNA set |
| popsize | Int | 50 | size of population |
| maxit | Int | 200 | iteration inner each time algorithm runs|
| optit | Int | 200 | number of optimize times |
| minval | Float | 2e-2 | min value to avoid zero, used in product|
| fitrev | Bool | false | whether used reversed H-M and Sm |
| norm | Bool|false | normalize objs before calculate distance between individuals and hyperplane|
| worstdef | Int | 0 | the method to choose the worst in DNA set, 0-product; 1-elucid distance|
| originpo | Bool | false | whether use original PO |
| cd | Bool | false | whether use crowding distance instead of knee point|

Notice the parameter `worstdef` has no real affect, 
code of elucid distance is removed cause it's too bad.

Please DO NOT use docker cause the docker file are not working yet.
