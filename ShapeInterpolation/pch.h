// Kullanmaya Başlama İpuçları: 
//   1. Dosyaları eklemek/yönetmek için Çözüm Gezgini penceresini kullanın
//   2. Kaynak denetimine bağlanmak için Takım Gezgini penceresini kullanın
//   3. Derleme çıktısını ve diğer iletileri görmek için Çıktı penceresini kullanın
//   4. Hataları görüntülemek için Hata Listesi penceresini kullanın
//   5. Yeni kod dosyaları oluşturmak için Projeye Git > Yeni Öğe ekle veya varolan kod dosyalarını projeye eklemek için Proje > Var Olan Öğeyi Ekle adımlarını izleyin
//   6. Bu projeyi daha sonra yeniden açmak için Dosya > Aç > Proje'ye gidip .sln uzantılı dosyayı seçin

#ifndef PCH_H
#define PCH_H

// TODO: add headers that you want to pre-compile here
#define _CRT_SECURE_NO_DEPRECATE
#define HAVE_SINGLEPRECISION_MATH //Eigen package does #include <complex> which conflicts with Inventor/C/basic.h definitions; so prevent them with this in files using include Inventor
#define _SECURE_SCL 0 //applications build on top of Open Inventor in release mode MUST define the _SECURE_SCL=0 preprocessor variable in order to avoid alignment problems. Using this mode reduce greatly memory footprint of STL objects and increase performance by a factor 2 to 10 depending on operations

#endif //PCH_H

/* ***************************** */
/*      MEMORY LEAK DETECTOR     */
/* ***************************** */
#define _CRTDBG_MAP_ALLOC
#include <cstdlib>
#include <crtdbg.h>

#ifdef _DEBUG
	#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
// Replace _NORMAL_BLOCK with _CLIENT_BLOCK if you want the
// allocations to be of _CLIENT_BLOCK type
#else
	#define DBG_NEW new
#endif
/* ***************************** */