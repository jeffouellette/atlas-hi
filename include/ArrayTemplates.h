#ifndef __ArrayTemplates_h__
#define __ArrayTemplates_h__

//namespace atlashi {

/**
 * Each of the following functions returns an n-dimensional array of T's with each entry set to 0.
 */
template <typename T>
T* Get1DArray (const int n1) {
  T* arr = new T[n1];
  for (int i = 0; i < n1; i++) {
   arr[i] = 0;
  }
  return arr;
}

template <typename T>
T** Get2DArray (const int n1, const int n2) {
  T** arr = new T*[n1];
  for (int i = 0; i < n1; i++) {
   arr[i] = Get1DArray <T> (n2);
  }
  return arr;
}

template <typename T>
T*** Get3DArray (const int n1, const int n2, const int n3) {
  T*** arr = new T**[n1];
  for (int i = 0; i < n1; i++) {
   arr[i] = Get2DArray <T> (n2, n3);
  }
  return arr;
}

template <typename T>
T**** Get4DArray (const int n1, const int n2, const int n3, const int n4) {
  T**** arr = new T***[n1];
  for (int i = 0; i < n1; i++) {
   arr[i] = Get3DArray <T> (n2, n3, n4);
  }
  return arr;
}

template <typename T>
T***** Get5DArray (const int n1, const int n2, const int n3, const int n4, const int n5) {
  T***** arr = new T****[n1];
  for (int i = 0; i < n1; i++) {
   arr[i] = Get4DArray <T> (n2, n3, n4, n5);
  }
  return arr;
}

template <typename T>
T****** Get6DArray (const int n1, const int n2, const int n3, const int n4, const int n5, const int n6) {
  T****** arr = new T*****[n1];
  for (int i = 0; i < n1; i++) {
   arr[i] = Get5DArray <T> (n2, n3, n4, n5, n6);
  }
  return arr;
}

/**
 * Each of the following functions deletes the given n-dimensional array of T's.
 * Since 2D, 3D, 4D,... arrays are always an array of pointers, they are all deleted the same way.
 * 1D arrays have their elements deleted only if the elements are pointers (e.g., to classes).
 */
template <typename T>
void Delete1DArray (T* &arr, const int n1, typename std::enable_if<std::is_pointer<T>::value, T>::type* = 0) {
  for (int i = 0; i < n1; i++) if (arr[i]) delete arr[i];
  delete[] arr;
  arr = NULL;
}

template <typename T>
void Delete1DArray (T* &arr, const int n1, typename std::enable_if<!std::is_pointer<T>::value, T>::type* = 0) {
  delete[] arr;
  arr = NULL;
}

template <typename T>
void Delete2DArray (T** &arr, const int n1, const int n2) {
  for (int i = 0; i < n1; i++) Delete1DArray (arr[i], n2);
  delete[] arr;
  arr = NULL;
}

template <typename T>
void Delete3DArray (T*** &arr, const int n1, const int n2, const int n3) {
  for (int i = 0; i < n1; i++) Delete2DArray (arr[i], n2, n3);
  delete[] arr;
  arr = NULL;
}

template <typename T>
void Delete4DArray (T**** &arr, const int n1, const int n2, const int n3, const int n4) {
  for (int i = 0; i < n1; i++) Delete3DArray (arr[i], n2, n3, n4);
  delete[] arr;
  arr = NULL;
}

template <typename T>
void Delete5DArray (T***** &arr, const int n1, const int n2, const int n3, const int n4, const int n5) {
   for (int i = 0; i < n1; i++) Delete4DArray (arr[i], n2, n3, n4, n5);
   delete[] arr;
   arr = NULL;
}

template <typename T>
void Delete6DArray (T****** &arr, const int n1, const int n2, const int n3, const int n4, const int n5, const int n6) {
   for (int i = 0; i < n1; i++) Delete5DArray (arr[i], n2, n3, n4, n5, n6);
   delete[] arr;
   arr = NULL;
}

//} // end namespace

#endif
