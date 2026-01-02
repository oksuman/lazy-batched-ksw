#pragma once

#include <variant>
#include "openfhe/pke/ciphertext.h"

namespace polycircuit
{

using Ciphertext = std::variant<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>;

class IComponent
{
public:
    virtual ~IComponent() = default;
    virtual Ciphertext evaluate() = 0;
    virtual Ciphertext evaluate_lazy() = 0;
};

} // namespace polycircuit
