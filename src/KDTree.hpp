// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include "Point.hpp"

template <typename value_type>
struct KDTreeNode
{
    value_type valor;
    KDTreeNode<value_type>* nodes[2];
    KDTreeNode(const value_type& value) {
        valor = value;
        nodes[0] = nodes[1] = 0;
    }
    KDTreeNode(value_type value, KDTreeNode<value_type>* left, KDTreeNode<value_type>* right)
    {
        valor = value;
        nodes[0] = left;
        nodes[1] = right;
    }
};

template <size_t N, typename ElemType>
class KDTree {
 public:
  typedef std::pair<Point<N>, ElemType> value_type;  //Importante------------------------------------------------

  KDTree();
  ~KDTree();

  KDTree(const KDTree &rhs);
  KDTree &operator=(const KDTree &rhs);

  size_t dimension() const;
  size_t size() const;
  bool empty() const;

  bool find(const Point<N>& pt, KDTreeNode<value_type>**& p) const;

  bool contains(const Point<N>& pt) const;
  void insert(const Point<N> &pt, const ElemType &value);
  ElemType &operator[](const Point<N> &pt);
  ElemType &at(const Point<N> &pt);
  const ElemType &at(const Point<N> &pt) const;

  void knn_aux(Point<N> key, KDTreeNode<value_type>* currentNode, KDTreeNode<value_type>*& guest, double& bestDistance, int nivel, std::vector<std::pair<ElemType, double> >& tipo) const;
  ElemType knn_value(const Point<N> &key, size_t k) const;
  std::vector<ElemType> knn_query(const Point<N> &key, size_t k) const;

 private:
  size_t dimension_;
  size_t size_;
  KDTreeNode<value_type>* root = nullptr;
};

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
    dimension_ = N;
    size_ = 0;
}

template <typename value_type>
void deleter(KDTreeNode<value_type>* node)
{
    if (node != nullptr) {
        deleter((node->nodes)[1]);
        deleter((node->nodes)[0]);
        delete node;
    }
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
    deleter(root);
    //root = nullptr;
}
template <typename value_type>
KDTreeNode<value_type>* copyNodes(const KDTreeNode<value_type>* node)
{
    if (node != nullptr){
        KDTreeNode<value_type>* nodeCopy = new KDTreeNode<value_type>(node->valor, copyNodes((node->nodes)[0]), copyNodes((node->nodes)[1]));
        return nodeCopy;
    }
    return nullptr;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree &rhs) {
    root = copyNodes(rhs.root);
    dimension_ = rhs.dimension_;
    size_ = rhs.size_;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType> &KDTree<N, ElemType>::operator=(const KDTree &rhs) {
    root = copyNodes(rhs.root);
    dimension_ = rhs.dimension_;
    size_ = rhs.size_;
    return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
  return N;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
  return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
  return !root;
}

template <size_t N, typename ElemType> 
bool KDTree<N, ElemType>::find(const Point<N>& pt, KDTreeNode<value_type>**& p) const {
    size_t i = 0;
    for (p = const_cast<KDTreeNode<value_type>**> (&root) ; *p && ((*p)->valor).first != pt; p = &((*p)->nodes[pt[i % N] > (((*p)->valor).first)[i % N]])) {//[adaptarse a cada sub nivel]
        i++;
    }
    return *p != 0;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N> &pt) const { 
    KDTreeNode<value_type>** p;
    if (!find(pt, p)) {
        return false;
    }
    return true;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N> &pt, const ElemType &value) {
    KDTreeNode<value_type>** p;
    if (!find(pt, p)) {
        value_type aux;
        aux.first = pt;
        aux.second = value;
        *p = new KDTreeNode<value_type>(aux);
        size_ += 1;
    }
    ((*p)->valor).second = value;
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::operator[](const Point<N> &pt) {
    KDTreeNode<value_type>** p;
    if (!find(pt, p)) {
        value_type aux;
        aux.first = pt;
        aux.second = size_;
        *p = new KDTreeNode<value_type>(aux);
        size_ += 1;
        return ((*p)->valor).second;
    }
    return ((*p)->valor).second;
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) {
    KDTreeNode<value_type>** p;
    if (find(pt, p)) {
        return ((*p)->valor).second;
    }
    throw std::out_of_range("");
}

template <size_t N, typename ElemType>
const ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) const {
    KDTreeNode<value_type>** p;
    if (find(pt, p)) {
        return ((*p)->valor).second;
    }
    throw std::out_of_range("");
}


template <size_t N, typename ElemType>
void KDTree<N, ElemType>::knn_aux(Point<N> key, KDTreeNode<value_type>* cNode, KDTreeNode<value_type>*& guest, double& optimo, int nivel, std::vector<std::pair<ElemType, double> >& tipo) const
{
    if (cNode == nullptr)
        return;
    double d = distance((cNode->valor).first, key);
    tipo.push_back(std::make_pair((cNode->valor).second, d));

    knn_aux(key, (cNode->nodes)[0], guest, optimo, ++nivel,tipo);
    knn_aux(key, (cNode->nodes)[1], guest, optimo, ++nivel,tipo);
}
template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N> &key, size_t k) const {
    if (k > size_)
        k = size_;
    std::vector<ElemType> values(k);
    values = knn_query(key, k);
    std::vector<ElemType> ElemV;
    ElemV.push_back(values[0]);
    for (size_t i = 1; i < k; i++) {
        bool ex = false;
        for (size_t j = 0; j < ElemV.size(); j++) {
            if (values[i] == ElemV[j])
                ex = true;
        }
        if(!ex)
            ElemV.push_back(values[i]);
    }
    std::vector<std::pair<ElemType, size_t> > cont;
    for (size_t i = 0; i < ElemV.size(); i++) {
        cont.push_back(std::make_pair(ElemV[i],0));
    }
    for (size_t i = 0; i < k; i++) {
        for (size_t j = 0; j < cont.size(); j++) {
            if (values[i] == cont[j].first)
                cont[j].second += 1;
        }
    }
    std::pair<ElemType, size_t> maxi = cont[0];
    for (size_t i = 1; i < cont.size(); i++) {
        if (maxi.second < cont[i].second)
            maxi = cont[i];
    }
    return maxi.first;
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N> &key, size_t k) const {
    std::vector<std::pair<ElemType, double> > tipo;
    KDTreeNode<value_type>* cNode = root;
    KDTreeNode<value_type>* guest = nullptr;
    double bestDistance = std::numeric_limits<double>::infinity();
    size_t nivel = 0;

    knn_aux(key, cNode, guest, bestDistance, nivel, tipo);
    std::sort(tipo.begin(), tipo.end(), [](const auto& x, const auto& y) { return x.second < y.second; });
    
    std::vector<ElemType> val;
    for (size_t i = 0; i < k; i++) {
        val.push_back(tipo[i].first);
    }
    return val;
}

#endif  // SRC_KDTREE_HPP_
