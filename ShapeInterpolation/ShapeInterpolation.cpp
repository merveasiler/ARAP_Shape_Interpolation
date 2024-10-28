// ShapeInterpolation.cpp : Bu dosya 'main' işlevi içeriyor. Program yürütme orada başlayıp biter.
//

/*
Yapilacaklar:
1) MeshTools.cpp'de 51 ve 52. satirlardaki crossProduct() ve computeLength() fonksiyonlari triangle'dan çıkarilip bagimsiz hâle getirilecek.
2) Mesh.h'daki Triangle struct'u icerisindeki crossProduct(), computeLength() ve normalize() fonksiyonları Triangle'ın icinden cikarilacak,
   ve herkese açık hale getirilecek. Bunları mesh.h'in icinde bagimsiz tanimladigimizda hata cikariyor (Mesh.h'i bir cok class icinde
   cagirdimizdan dolayi). Ayrica, bu fonksiyonları tüm fonksiyonlarin icinde tekrar tekrar imlement etmektense bir kez yazilmis olanlari
   cagirilacak.
3) libigl duzenlenecek. Add-link etme islemleri yani.
7) Sparse matrix'ler symmetric formda olduklari icin islem sayisini yariya dusurebilirsin.
8) Float islemleri kusurat hassasiyetinde uydurmalar var. Hassasiyete limit getirebilirsin.
9) SVD garip calisiyor. Rotation kismi orthogonal cikmiyor. Manuel SVD yaptirabilirsin.
*/

//ShapeInterpolation.cpp: Konsol uygulamasının giriş noktasını tanımlar.
//

#include "pch.h"
#include "Transformation.h"
#include "ShapeInterpolation.h"

int main(int argc, char * argv[])
{
	
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_DEBUG);
	{
		// INTERPOLATE:
		// Example: interpolate DB_Hand/hand_gestures_off/0.off DB_Hand/hand_gestures_off/10.off 0.5 DB_Hand/method_alexa/0to10at0.5.off
		if (argv[1][0] == 'i')
			interpolate(argv[2], argv[3], stof(argv[4], 0), argv[5]);

		// DRAW:
		// Example:	draw DB_Kids/method_alexa/11to15at0.5.off
		else if (argv[1][0] == 'd')
			drawMeshToScene(argv[2]);

		// UNDEFINED:
		else
			cout << "Undefined Operation!" << endl;
	}
	_CrtDumpMemoryLeaks();

	return 0;

}

